import cv2
import yaml
import fire
import os
import csv
import numpy as np
from math import dist, sqrt, ceil, atan2
from typing import Tuple, List


class MapManager:
    def __init__(self, map_yaml_file) -> None:
        self.map_yaml_file = map_yaml_file
        self.window_name = "Course Image"
        self.map_img, self.map_yaml_obj = self.read_map(map_yaml_file)
        self.map_height, self.map_width = self.map_img.shape[:2]

    @staticmethod
    def read_map(map_yaml_file) -> Tuple[np.ndarray, dict]:
        yaml_obj = get_yaml_obj(map_yaml_file)
        map_dir_path = os.path.dirname(map_yaml_file)
        map_map_img_path = os.path.join(map_dir_path, yaml_obj["image"])
        map_img = cv2.imread(map_map_img_path, cv2.IMREAD_COLOR)
        return map_img, yaml_obj

    def pixel2metric_point(self, p_point) -> list:
        m_point = [
            self.map_yaml_obj["origin"][0]
            + self.map_yaml_obj["resolution"] / 2
            + p_point[0] * self.map_yaml_obj["resolution"],
            self.map_yaml_obj["origin"][1]
            + self.map_yaml_obj["resolution"] / 2
            + (self.map_height - p_point[1] - 1) * self.map_yaml_obj["resolution"]
        ]
        return m_point

    def pixel2metric_line(self, p_line) -> list:
        m_line = []
        for p_point in p_line:
            m_point = self.pixel2metric_point(p_point)
            m_line.append(m_point)
        return m_line

    def pixel2metric_lines(self, p_lines) -> np.ndarray:
        m_lines = []
        for p_line in p_lines:
            m_line = self.pixel2metric_line(p_line)
            m_lines.append(m_line)
        return np.array(m_lines)

    def metric2pixel_point(self, m_point) -> list:
        p_point = [
            (m_point[0] - self.map_yaml_obj["origin"][0] - self.map_yaml_obj["resolution"] / 2) / self.map_yaml_obj["resolution"],
            (self.map_height - 1 - (m_point[1] - self.map_yaml_obj["origin"][1] - self.map_yaml_obj["resolution"] / 2) / self.map_yaml_obj["resolution"])
        ]
        return p_point

    def metric2pixel_line(self, m_line) -> list:
        p_line = []
        for m_point in m_line:
            p_point = self.metric2pixel_point(m_point)
            p_line.append(p_point)
        return p_line

    def metric2pixel_lines(self, m_lines) -> np.ndarray:
        p_lines = []
        for m_line in m_lines:
            p_line = self.metric2pixel_line(m_line)
            p_lines.append(p_line)
        return np.array(p_lines)

    def metric2pixel_length(self, m_length) -> float:
        p_length = m_length / self.map_yaml_obj["resolution"]
        return p_length

    def resize_image(self, img) -> np.ndarray:
        img_height = img.shape[0]
        img_width = img.shape[1]
        if img_height > 1080:
            scale = 1080 / img_height
            img = cv2.resize(img, dsize=None, fx=scale, fy=scale)
        if img_width > 1920:
            scale = 1920 / img_width
            img = cv2.resize(img, dsize=None, fx=scale, fy=scale)
        return img

    def show_image(self) -> None:
        resized_map_img = self.resize_image(self.map_img)
        cv2.imshow(self.window_name, resized_map_img)
        try:
            while True:
                if (
                    cv2.waitKey(1) & 0xFF == ord("q")
                    or cv2.getWindowProperty(self.window_name, cv2.WND_PROP_VISIBLE)
                    <= 0
                ):
                    break
        except KeyboardInterrupt:
            cv2.destroyAllWindows()


class ContourDetector(MapManager):
    def __init__(self, map_yaml_file) -> None:
        super().__init__(map_yaml_file)

    def get_angle(self, vector1, vector2) -> float:
        angle = atan2(vector1[1], vector1[0]) - atan2(vector2[1], vector2[0])
        if angle <= -np.pi:
            angle += 2 * np.pi
        elif angle > np.pi:
            angle -= 2 * np.pi
        return angle

    def get_course_contours(self, metric_crop_length) -> list:
        p_crop_length = self.metric2pixel_length(metric_crop_length)
        map_img_gray = cv2.cvtColor(self.map_img, cv2.COLOR_BGR2GRAY)  # convert to gray scale
        map_img_bin = cv2.threshold(map_img_gray, 250, 255, cv2.THRESH_BINARY)[1]  # biniarize
        map_img_eroded = cv2.erode(map_img_bin, cv2.getStructuringElement(cv2.MORPH_RECT, (3,3)), iterations=int(p_crop_length)) # erode
        contours, hierarchy = cv2.findContours(map_img_eroded, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_TC89_L1)  # not_detected contours

        # add contour area data to hierarchy
        hierarchy = hierarchy.reshape(len(hierarchy[0]), 4)
        hierarchy_w_area = []
        for idx, (contour, data) in enumerate(zip(contours, hierarchy)):
            data_w_area = [idx, cv2.contourArea(contour), data]
            hierarchy_w_area.append(data_w_area)

        # sort contours by area
        hierarchy_w_area.sort(key=lambda x: x[1], reverse=True)

        outer_idx, inner_idx = None, None
        for (idx, area, data) in hierarchy_w_area:
            if (cv2.contourArea(contour) < (self.map_height-p_crop_length-2)*(self.map_width-p_crop_length-2) # remove image edge contour
                and data[2] >= 0): # remove contours with no pairs
                outer_idx = idx
                inner_idx = data[2]

        # get outer and inner contours of the course
        pixel_contours = [
            contours[outer_idx],  # outer contour
            contours[inner_idx],  # inner contour
        ]
        pixel_contours = [
            contour.reshape(len(contour), 2) for contour in pixel_contours # reshape to points array
        ]

        # convert pixel coordinates to metric coordinates
        metric_contours = []
        for contour in pixel_contours:
            metric_contour = self.pixel2metric_line(contour)
            metric_contours.append(np.array(metric_contour))

        return metric_contours

    def draw_contours(self, img, contour_list, bgr=(0, 0, 255), thickness=1) -> np.ndarray:
        m_contours = []
        for contour in contour_list:
            m_contour = self.metric2pixel_line(contour)
            m_contour = np.array(m_contour, dtype=np.int32)
            m_contours.append(m_contour)
        img = cv2.drawContours(img, m_contours, -1, bgr, thickness)
        return img

    def draw_points(self, img, point_list, bgr=(0, 255, 0), size=2, thickness=1) -> np.ndarray:
        point_list = self.metric2pixel_line(point_list)
        for point in point_list:
            cv2.drawMarker(
                img,
                [int(point[0]), int(point[1])],
                bgr,
                cv2.MARKER_CROSS,
                size,
                thickness,
            )
        return img


class LateralLineCreator(ContourDetector):
    def __init__(
        self,
        map_yaml_file
    ) -> None:
        super().__init__(map_yaml_file)

    def get_nearest_index(
        self, ref_point, points_list, prev_nearest_idx=None, search_length=None
    ) -> int:
        """Returns the index of the nearest point to the reference point in the given point list."""
        min_idx = 0
        if (prev_nearest_idx is not None and search_length is not None): # if search range is specified
            start_idx = (prev_nearest_idx - search_length) % len(points_list)
            end_idx = (prev_nearest_idx + search_length) % len(points_list)
            if start_idx > end_idx:
                start_idx -= len(points_list)
            clamp_list = []
            for i in range(start_idx, end_idx + 1):
                clamp_list.append(points_list[i])
            points_list = clamp_list
            min_idx += start_idx
        min_idx += min(enumerate(points_list), key=lambda x: dist(ref_point, x[1]))[0]
        return min_idx

    def get_point2line_position(self, point, line) -> int:
        """
        Tells whether the nearest point on the line to the given point is on the end or in the middle of the line.
        Returns: -1 if the nearest point is line[0], 1 if the nearest point is line[1], 0 if the nearest point is in the middle of the line.
        """
        dot_start = (point[0] - line[0][0]) * (line[1][0] - line[0][0]) + (
            point[1] - line[0][1]
        ) * (line[1][1] - line[0][1])
        dot_end = (point[0] - line[1][0]) * (line[0][0] - line[1][0]) + (
            point[1] - line[1][1]
        ) * (line[0][1] - line[1][1])
        if dot_start <= 0:
            return -1
        elif dot_end <= 0:
            return 1
        else:
            return 0

    def get_point2line_distance(self, point, line) -> float:
        """
        Returns the distance between the given point and the nearest point on the line.
        Values are calculated assuming an infinite length line.
        """
        return abs(
            (line[1][0] - line[0][0]) * (line[0][1] - point[1])
            - (line[0][0] - point[0]) * (line[1][1] - line[0][1])
        ) / dist(line[0], line[1])

    def get_line_data_list(self, ref_point, line_list) -> Tuple[List[float], List[float]]:
        dist_list = [] # distance from the reference point to the nearest point on the line
        proj_list = [] # projection ratio of the vector to the reference point on the line
        for idx, line in enumerate(line_list):
            pos = self.get_point2line_position(ref_point, line)
            if pos == 0:
                line_dist = self.get_point2line_distance(ref_point, line)
                dist_list.append(line_dist)
                ref_dist = dist(line[0], ref_point)
                if ((ref_dist * ref_dist - line_dist * line_dist) < 0):
                    proj = 0
                else:
                    proj = sqrt(ref_dist * ref_dist - line_dist * line_dist) / dist(line[0], line[1])
                proj_list.append(proj)
            elif pos == -1:
                dist_list.append(dist(line[0], ref_point))
                proj_list.append(0)
            elif pos == 1:
                dist_list.append(dist(line[1], ref_point))
                proj_list.append(1)
        return dist_list, proj_list

    def get_lateral_lines(self, contour_list) -> list:
        """
        Returns a list of lateral lines of the race track.
        """
        if len(contour_list) != 2:
            raise ValueError(
                f"Invalid numbers of contours: {len(contour_list)}. (Expected: 2)"
            )

        lateral_line_list = []
        for i in range(2):  # outer and inner contours
            ref_contour = contour_list[i-1]
            srch_contour = contour_list[i]
            for ref_idx, ref_point in enumerate(ref_contour):
                line_list = []
                for idx in range(len(srch_contour)):
                    line_list.append([
                        srch_contour[idx-1],
                        srch_contour[idx],
                    ])
                dist_list, proj_list = self.get_line_data_list(ref_point, line_list) # get distances and projection ratios
                line_idx, (line, line_dist, proj) = min(
                    enumerate(zip(line_list, dist_list, proj_list)),
                    key=lambda x: x[1][1],
                )  # get nearest line
                nearest_point = (
                    line[0] + (line[1] - line[0]) * proj
                )  # get nearest point on the line
                if i == 0: # inner contour is the reference
                    lateral_line_list.append([nearest_point, ref_point])
                else: # outer contour is the reference
                    lateral_line_list.append([ref_point, nearest_point])
        return lateral_line_list

    def center_point(self, line) -> np.ndarray:
        """
        Returns the center point of the given line.
        """
        return np.array([
                (line[0][0] + line[1][0]) / 2,
                (line[0][1] + line[1][1]) / 2
        ])

    def sort_lines(self, line_list, contour_list, is_clockwise=True):
        # get contour lines
        contour_line_list = []
        for idx in range(len(contour_list[0])):
            contour_line_list.append([
                contour_list[0][idx-1],
                contour_list[0][idx]
            ])

        # prepare list for sorting
        extended_list = []
        for line in line_list:
            extended_list.append([line, self.center_point(line), 0])

        # remove lines that are not on the contour
        for idx in reversed(range(len(extended_list))):
            not_detected = True
            for contour_line in contour_line_list:
                not_detected = not_detected and not self.detect_overlaps(extended_list[idx][0], contour_line)
            if not_detected:
                extended_list.pop(idx)

        # get first line
        sorted_list = []
        overlap_list = []
        for i, contour_line in enumerate(contour_line_list):
            for j in reversed(range(len(extended_list))):
                if self.detect_overlaps(extended_list[j][0], contour_line):
                    overlap_list.append(extended_list.pop(j))
            if len(overlap_list) == 1:
                sorted_list.append(overlap_list[0])
            elif len(overlap_list) > 1:
                for k, data in enumerate(overlap_list):
                    overlap_list[k][2] = dist(contour_line_list[i+1][0], data[1])
                overlap_list.sort(key=lambda x: x[2], reverse=True)
                while len(overlap_list) > 0:
                    sorted_list.append(overlap_list.pop(0))
            if len(sorted_list) > 0:
                break

        # add the rest of the lines in order
        for contour_line in contour_line_list:
            overlap_list = []
            for idx in reversed(range(len(extended_list))):
                if self.detect_overlaps(extended_list[idx][0], contour_line):
                    overlap_list.append(extended_list.pop(idx))
            if len(overlap_list) == 1:
                sorted_list.append(overlap_list[0])
            elif len(overlap_list) > 1:
                for idx, data in enumerate(overlap_list):
                    overlap_list[idx][2] = dist(sorted_list[-1][1], data[1])
                overlap_list.sort(key=lambda x: x[2])
                while len(overlap_list) > 0:
                    sorted_list.append(overlap_list.pop(0))

        # get simple list of lines
        sorted_line_list = [data[0] for data in sorted_list]

        # check direction
        angle_sum = 0
        for idx, line in enumerate(sorted_line_list):
            vector1 = sorted_line_list[idx][1] - sorted_line_list[idx-1][1]
            vector2 = sorted_line_list[idx-1][1] - sorted_line_list[idx-2][1]
            angle_sum += self.get_angle(vector1, vector2)
        if is_clockwise and angle_sum > 0:
            sorted_line_list = sorted_line_list[::-1]
        elif not is_clockwise and angle_sum < 0:
            sorted_line_list = sorted_line_list[::-1]

        return sorted_line_list

    def detect_overlaps(self, line1, line2) -> bool:
        """
        Returns True if the given lines overlap.
        This process considers the endpoints.
        """
        ax = line1[0][0]
        ay = line1[0][1]
        bx = line1[1][0]
        by = line1[1][1]
        cx = line2[0][0]
        cy = line2[0][1]
        dx = line2[1][0]
        dy = line2[1][1]
        s = (ax - bx) * (cy - ay) - (ay - by) * (cx - ax)
        t = (ax - bx) * (dy - ay) - (ay - by) * (dx - ax)
        if s*t > 0:
            return False
        elif s*t == 0:
            # one edge on line
            if t != 0:
                if (
                    (cx < ax and cx < bx) or (cx > ax and cx > bx)
                    or (cy < ay and cy < by) or (cy > ay and cy > by)
                ):
                    return False
            # one edge on line
            if s != 0:
                if (
                    (dx < ax and dx < bx) or (dx > ax and dx > bx)
                    or (dy < ay and dy < by) or (dy > ay and dy > by)
                ):
                    return False
            # both edges on line
            if (
                (cx < ax and cx < bx and dx < ax and dx < bx)
                or (cx > ax and cx > bx and dx > ax and dx > bx)
                or (cy < ay and cy < by and dy < ay and dy < by)
                or (cy > ay and cy > by and dy > ay and dy > by)
            ):
                return False

        s = (cx - dx) * (ay - cy) - (cy - dy) * (ax - cx)
        t = (cx - dx) * (by - cy) - (cy - dy) * (bx - cx)
        if s*t > 0:
            return False
        elif s*t == 0:
            # one edge on line
            if t != 0:
                if (
                    (ax < cx and ax < dx) or (ax > cx and ax > dx)
                    or (ay < cy and ay < dy) or (ay > cy and ay > dy)
                ):
                    return False
            # one edge on line
            if s != 0:
                if (
                    (bx < cx and bx < dx) or (bx > cx and bx > dx)
                    or (by < cy and by < dy) or (by > cy and by > dy)
                ):
                    return False
            # both edges on line
            if (
                (ax < cx and ax < dx and bx < cx and bx < dx)
                or (ax > cx and ax > dx and bx > cx and bx > dx)
                or (ay < cy and ay < dy and by < cy and by < dy)
                or (ay > cy and ay > dy and by > cy and by > dy)
            ):
                return False
        return True

    def detect_mid_overlaps(self, line1, line2) -> bool:
        """
        Returns True if the given lines overlap.
        This process does not consider the endpoints.
        """
        ax = line1[0][0]
        ay = line1[0][1]
        bx = line1[1][0]
        by = line1[1][1]
        cx = line2[0][0]
        cy = line2[0][1]
        dx = line2[1][0]
        dy = line2[1][1]
        s = (ax - bx) * (cy - ay) - (ay - by) * (cx - ax)
        t = (ax - bx) * (dy - ay) - (ay - by) * (dx - ax)
        if s*t > 0:
            return False
        elif s*t == 0:
            # one edge on line
            if s != 0 or t != 0:
                return False
            # both edges on line
            if (
                (cx <= ax and cx <= bx and dx <= ax and dx <= bx)
                or (cx >= ax and cx >= bx and dx >= ax and dx >= bx)
            ):
                return False

        s = (cx - dx) * (ay - cy) - (cy - dy) * (ax - cx)
        t = (cx - dx) * (by - cy) - (cy - dy) * (bx - cx)
        if s*t > 0:
            return False
        elif s*t == 0:
            # one edge on line
            if s != 0 or t != 0:
                return False
            # both edges on line
            if (
                (ax <= cx and ax <= dx and bx <= cx and bx <= dx)
                or (ax >= cx and ax >= dx and bx >= cx and bx >= dx)
                or (ay <= cy and ay <= dy and by <= cy and by <= dy)
                or (ay >= cy and ay >= dy and by >= cy and by >= dy)
            ):
                return False
        return True

    def remove_self_overlaps(self, line_list: list) -> np.ndarray:
        overlap_array : np.ndarray = np.zeros((len(line_list),len(line_list)), dtype=np.uint8)
        for i in range(len(line_list)):
            for j in range(i + 1, len(line_list)):
                if self.detect_mid_overlaps(line_list[i], line_list[j]):
                    overlap_array[i][j] += 1
        overlap_array = overlap_array + overlap_array.T
        mask = np.full(len(line_list), True, dtype=bool)
        try:
            while True:
                overlap_count = sum(overlap_array)
                if sum(overlap_count) == 0:
                    break
                max_idx = np.argsort(overlap_count)[-1]
                mask[max_idx] = False
                overlap_array[:, max_idx] = 0
                overlap_array[max_idx, :] = 0
        except KeyboardInterrupt:
            pass
        line_list = np.array(line_list)[mask]
        return line_list

    def remove_contour_overlaps(self, line_list, contour_list) -> list:
        """
        Returns a list of lines with overlapping lines with contours removed.
        """
        filtered_line_list = []
        for line in line_list:
            delete_flag = False
            for contour in contour_list:
                overlap_count = 0
                for idx, point in enumerate(contour):
                    contour_line = [contour[idx], contour[idx-1]]
                    if self.detect_overlaps(line, contour_line):
                        overlap_count += 1
                        if overlap_count >= 3:
                            delete_flag = True
                            break
                if delete_flag:
                    break
            if not delete_flag:
                filtered_line_list.append(line)
        return filtered_line_list

    def fill_gaps(self, line_list, max_gap) -> list:
        """
        Returns a list of lines with gaps filled.
        Unit of the gap is the same as the given line list.
        """
        filled_line_list = []
        for idx, line in enumerate(line_list):
            filled_line_list.append(line_list[idx - 1])
            start_dist = dist(line[0], line_list[idx - 1][0])
            end_dist = dist(line[1], line_list[idx - 1][1])
            if start_dist > max_gap or end_dist > max_gap:
                gap_ratio = max(ceil(start_dist / max_gap), ceil(end_dist / max_gap))
                for i in range(1, gap_ratio):
                    filled_line_list.append(
                        [
                            line_list[idx - 1][0]
                            + (line[0] - line_list[idx - 1][0]) * i / gap_ratio,
                            line_list[idx - 1][1]
                            + (line[1] - line_list[idx - 1][1]) * i / gap_ratio,
                        ]
                    )
        return filled_line_list

    def thin_out(self, line_list, min_gap) -> list:
        """
        Returns a list of lines thinned out based on minimum gap.
        Unit of the gap is the same as the given line list.
        """
        thinned_line_list = [line_list[0]]
        for idx, line in enumerate(line_list):
            d_1 = self.get_point2line_distance(line[0], thinned_line_list[-1])
            d_2 = self.get_point2line_distance(line[1], thinned_line_list[-1])
            if ((d_1 >= min_gap or d_2 >= min_gap) and d_1 > 0 and d_2 > 0):
                thinned_line_list.append(line)
        return thinned_line_list

    def crop_edges(self, line_list, crop_length) -> list:
        """
        Returns a list of lines cropped at both ends.
        Unit of the crop length is the same as the given line list.
        """
        cropped_line_list = []
        for line in line_list:
            cropped_line_list.append([
                    line[0] + (line[1] - line[0]) * crop_length / dist(line[0], line[1]),
                    line[1] + (line[0] - line[1]) * crop_length / dist(line[0], line[1]),
                ]
            )
        return cropped_line_list

    def draw_line(self, img, line, bgr=(0, 0, 255), thickness=1) -> np.ndarray:
        line = self.metric2pixel_line(line)
        cv2.line(
            img,
            [int(line[0][0]), int(line[0][1])],
            [int(line[1][0]), int(line[1][1])],
            bgr,
            thickness,
        )
        return img

    def draw_lines(self, img, line_list, bgr=(0, 0, 255), thickness=1) -> np.ndarray:
        line_list = self.metric2pixel_lines(line_list)
        for line in line_list:
            cv2.line(
                img,
                [int(line[0][0]), int(line[0][1])],
                [int(line[1][0]), int(line[1][1])],
                bgr,
                thickness,
            )
        return img

    def draw_polyline(self, img, polyline, bgr=(0, 0, 255), thickness=1) -> np.ndarray:
        polyline = self.metric2pixel_line(polyline)
        cv2.polylines(img, [np.array(polyline).astype(np.int32)], True, bgr, thickness)
        return img


def get_yaml_obj(yaml_path) -> dict:
    with open(yaml_path, "r") as stream:
        try:
            yaml_obj = yaml.safe_load(stream)
        except yaml.YAMLError as e:
            print(e)
    return yaml_obj


def check_paths(map_name: str) -> Tuple[str, str, str]:
    cdp = os.path.dirname(os.path.abspath(__file__))

    map_dir = os.path.join(cdp, "map")
    if not os.path.isdir(map_dir):
        raise FileNotFoundError(f"Map directory not found: {map_dir}")
    else:
        map_yaml_file = os.path.join(map_dir, map_name + ".yaml")    
        if not os.path.isfile(map_yaml_file):
            raise FileNotFoundError(f"Map yaml file not found: {map_yaml_file}")

    data_dir = os.path.join(cdp, "data")
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    out_dir = os.path.join(data_dir, map_name)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    config_dir = os.path.join(cdp, "config")
    if not os.path.isdir(config_dir):
        raise FileNotFoundError(f"Config directory not found: {config_dir}")
    config_file = os.path.join(config_dir, "path_optimizer.yaml")
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")

    return map_yaml_file, out_dir, config_file


def detect_course(map_name: str = None) -> Tuple[LateralLineCreator, List[List[List[float]]]]:
    if map_name is None:
        raise ValueError("map_name is None. \nPlease specify map_name.\n\n\tpython3 course_detection.py [map_name]\n\tpython3 path_optimizer.py [map_name]\n")

    # check if map and data directories exist
    map_yaml_file, out_dir, config_file = check_paths(map_name)

    # load course prep parameters
    path_yaml_obj = get_yaml_obj(config_file)
    is_clockwise = path_yaml_obj["course_prep"]["is_clockwise"]
    m_max_gap = path_yaml_obj["course_prep"]["interpolation_max_gap"]
    m_min_gap = path_yaml_obj["course_prep"]["thin_out_min_gap"]
    m_crop_length = path_yaml_obj["course_prep"]["crop_length"]
    bgr_list = path_yaml_obj["visualization"]["bgr_list"]

    # preperate course
    lat_line_creator = LateralLineCreator(
        map_yaml_file
    )

    print("Detecting course...")
    contour_list = lat_line_creator.get_course_contours(m_crop_length)
    print("Done.")
    print("Generating lateral lines...")
    line_list_raw = lat_line_creator.get_lateral_lines(contour_list)
    line_list_no_self_overlap = lat_line_creator.remove_self_overlaps(line_list_raw)
    line_list_no_overlap = lat_line_creator.remove_contour_overlaps(line_list_no_self_overlap, contour_list)
    line_list_sorted = lat_line_creator.sort_lines(line_list_no_overlap, contour_list, is_clockwise)
    line_list_dense = lat_line_creator.fill_gaps(line_list_sorted, m_max_gap)
    line_list_cropped = lat_line_creator.crop_edges(line_list_dense, 1e-2)
    line_list_thinned = lat_line_creator.thin_out(line_list_cropped, m_min_gap)
    center_points = [
        [(line[0][0] + line[1][0]) / 2, (line[0][1] + line[1][1]) / 2]
        for line in line_list_dense
    ]

    map_img = lat_line_creator.map_img.copy()
    map_img_w_contour_points = map_img.copy()
    for contour in contour_list:
        map_img_w_contour_points = lat_line_creator.draw_points(map_img_w_contour_points, contour)

    # save image with contours
    img = lat_line_creator.draw_contours(map_img_w_contour_points.copy(), contour_list)
    file_path = os.path.join(out_dir, f"{map_name}_w_contours.png")
    cv2.imwrite(file_path, img)

    # save image with lateral lines
    img = lat_line_creator.draw_lines(map_img_w_contour_points.copy(), line_list_raw)
    file_path = os.path.join(out_dir, f"{map_name}_w_lateral_lines_raw.png")
    cv2.imwrite(file_path, img)

    # save image with overlap free lateral lines
    img = lat_line_creator.draw_lines(map_img_w_contour_points.copy(), line_list_no_overlap)
    file_path = os.path.join(out_dir, f"{map_name}_w_overlap_free_lines.png")
    cv2.imwrite(file_path, img)

    # save image with adjusted lateral lines
    img = lat_line_creator.draw_lines(map_img_w_contour_points.copy(), line_list_thinned)
    file_path = os.path.join(out_dir, f"{map_name}_w_lateral_lines_adjusted.png")
    cv2.imwrite(file_path, img)

    # save image with center lines
    img = map_img.copy()
    img = lat_line_creator.draw_polyline(map_img.copy(), center_points)
    file_path = os.path.join(out_dir, f"{map_name}_w_center_line.png")
    cv2.imwrite(file_path, img)

    # save course data as csv
    file_path = os.path.join(out_dir, f"{map_name}_course_info.csv")
    with open(file_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["# x_m", "y_m", "w_tr_right", "w_tr_left", "outer_x", "outer_y", "inner_x", "inner_y"])
        for center_point, lat_line in zip(center_points, line_list_dense):
            width = dist(lat_line[0], lat_line[1])
            writer.writerow([center_point[0], center_point[1], width/2,  width/2, lat_line[0][0], lat_line[0][1], lat_line[1][0], lat_line[1][1]])

    print("Done.")

    return lat_line_creator, line_list_thinned


def main(map_name: str = None) -> Tuple[LateralLineCreator, List[List[List[float]]]]:
    detect_course(map_name)


if __name__ == "__main__":
    fire.Fire(main)
