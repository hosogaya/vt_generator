#include <vt_generator/csv/csv_reader_opt_curvature.hpp>
#include <vt_generator/csv/csv_writer.hpp>
#include <vt_generator/type.hpp>
#include <vt_generator/util.hpp>

using namespace vt_generator;
int main()
{
    std::string file_name = "../data/exp20240422/opt_path_curvature.csv";
    csv::Reader reader(file_name);

    // cal center pos
    std::vector<Vector2> center;
    for (int i=0; i<reader.size(); ++i)
    {
        center.push_back(Vector2{reader.xm()[i], reader.ym()[i]});
    }

    // cal normal vectors
    std::vector<Vector2> normal;
    for (int i=0; i<reader.size(); ++i)
    {
        normal.push_back(calNormalVector2(reader.xm(), reader.ym(), i)); // towards left hand
    }

    std::vector<Vector2> outer_line;
    std::vector<Vector2> inner_line;
    std::vector<Scalar> outer_w, inner_w;
    for (int i=0; i<reader.size(); ++i)
    {
        outer_line.push_back(getClosestPointOnLine(-normal[i], center[i], reader.outer_x(), reader.outer_y())); // right
        inner_line.push_back(getClosestPointOnLine(normal[i], center[i], reader.inner_x(), reader.inner_y())); // left
        outer_w.push_back((outer_line.at(i) - center.at(i)).norm());
        inner_w.push_back((inner_line.at(i) - center.at(i)).norm());
    }

    
    std::vector<Scalar> outer_x, outer_y, inner_x, inner_y;
    for (int i=0; i<reader.size(); ++i)
    {
        // outer_x.push_back(reader.center_x()[i] - normal[i].x());
        // outer_y.push_back(reader.center_y()[i] - normal[i].y());
        outer_x.push_back(outer_line.at(i).x());
        outer_y.push_back(outer_line.at(i).y());
        inner_x.push_back(inner_line.at(i).x());
        inner_y.push_back(inner_line.at(i).y());
    }
    csv::Writer writer("line_modified.csv");
    writer.writeResult(reader.xm(), reader.ym(), outer_w, inner_w, 
                      reader.center_x(), reader.center_y(), outer_x, outer_y, inner_x, inner_y, 
                    reader.curvature(), reader.ref_v(), reader.xm().size());


    return 1;
}   