//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"

substrate::substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces, transform &t, std::function<double(double)> strain_f) : _tree()
{
    _transform = std::move(t);
    _strain = strain_f;
    _myocytes.reserve(myo_vertices.size());
    _points.reserve(myo_vertices.size());
    int i=0;
    for (auto &myo_vertices_i : myo_vertices)
    {
        polyhedronSet myo(myo_vertices_i, myo_faces[i]);
        Eigen::Vector3d centroid = myo.getUnstrainedCentroid();
        Kernel::Point_3 centroid_cgal(centroid(0), centroid(1), centroid(2));
        _points.emplace_back(centroid_cgal);
        _map_centroid_to_polygon[centroid_cgal] = i; 
        _myocytes.push_back(std::move(myo));
        i++;
    }
    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
}

substrate::substrate(const std::vector<std::string>& polygons, transform &t, std::function<double(double)> strain_f){

    _transform = std::move(t);
    _strain = strain_f;

    //iterate over polygons
    int i=0;
    auto start = std::chrono::high_resolution_clock::now(); // start timer

    for (const auto &poly: polygons){
        polyhedronSet myo(poly);
        Eigen::Vector3d centroid = myo.getUnstrainedCentroid();
        Kernel::Point_3 centroid_cgal(centroid(0), centroid(1), centroid(2));
        _points.emplace_back(centroid_cgal);
        _map_centroid_to_polygon[centroid_cgal] = i;
        _myocytes.push_back(std::move(myo));       
        i++;
    }
    
    Eigen::Vector3d centroid_myo = _myocytes[0].getUnstrainedCentroid();

    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
    auto end = std::chrono::high_resolution_clock::now(); // end timer
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time set substrate: " << elapsed.count() << " seconds" << std::endl;
}

void substrate::write_ply_surface_mesh(const Mesh &mesh, const std::string &filename)
{

    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << mesh.number_of_vertices() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << mesh.number_of_faces() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &v : mesh.vertices())
    {
        const auto &point = mesh.point(v);
        out << point.x() << " " << point.y() << " " << point.z() << "\n";
    }

    // Write faces
    for (const auto &f : mesh.faces())
    {
        out << "3";
        CGAL::Vertex_around_face_circulator<Mesh> vcirc(mesh.halfedge(f), mesh), done(vcirc);
        do
        {
            out << " " << (*vcirc).idx();
        } while (++vcirc != done);
        out << "\n";
    }

    out.close();
}

void substrate::write_ply_polyhedron(polygon &poly, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    // Write header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << poly.getVertices().size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << poly.getFaces().size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (auto &vertex : poly.getVertices())
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    for (auto &face : poly.getFaces())
    {
        out << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
    }
    out.close();
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>> substrate::intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const
{
    // We will place the point in the center of the block as the tree is built with the centroid of the polygons

    // TODO: This will have to be modified when considering a full geometry. Maybe input into to the Ktree the 8 points of bounding box?
    // TODO: CHECK THAT THIS WORKS CURRENTLY, I DO NOT SEE WHY IT SHOULDNT!
    Kernel::Point_3 query((double)point(0), (double)point(1), (double)point(2));

    Neighbor_search search(*_tree, query, 5);

    double min_distance_children = std::numeric_limits<double>::max();
    int index_polygon = -1;
    Eigen::Vector3d normal_vector = Eigen::Vector3d::Zero();
    
    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        //Pass in total_time here
    
        int index_polygon = _map_centroid_to_polygon.at(it->first);
        const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
        boost::optional<std::pair<int, double>> intersection = curr_poly.intersection(point, step);
        if (intersection)
        {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            Kernel::FT distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;

            if (distance_to_face < min_distance_children)
            {
                min_distance_children = CGAL::to_double(distance_to_face);
                index_polygon = _map_centroid_to_polygon.at(it->first);
                normal_vector = curr_poly.getNormalVector(index_face);
            }
        }
    }

    if (min_distance_children == std::numeric_limits<double>::max())
    {
        return boost::optional<std::tuple<int, double, Eigen::Vector3d>>();
    }
    
    return std::make_tuple(index_polygon, min_distance_children, normal_vector);
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>>  substrate::intersectionBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const
{
    
    // TODO: If we want reflecting boundaries in the block we need to output the normal vector of the face intersected
    CGAL::Bbox_3 solid_box = _transform.getBlockCurrentTime(index_sequence).getSolidBbox();
    CGAL::Point_3<Kernel> p1(point(0), point(1), point(2));
    CGAL::Point_3<Kernel> p2(point(0) + step(0), point(1) + step(1), point(2) + step(2));
    CGAL::Segment_3<Kernel> segment(p1, p2);

    //Check if p1 point is outside of the solid_box
    if (p1.x() < solid_box.xmin() || p1.x() > solid_box.xmax() ||
        p1.y() < solid_box.ymin() || p1.y() > solid_box.ymax() ||
        p1.z() < solid_box.zmin() || p1.z() > solid_box.zmax())
    {
        throw std::runtime_error("how did we get here? p1 is outside the solid_box");
    }

    // Check if point p2 is inside the solid_box without the boundaries
    if (p2.x() > solid_box.xmin() && p2.x() < solid_box.xmax() &&
        p2.y() > solid_box.ymin() && p2.y() < solid_box.ymax() &&
        p2.z() > solid_box.zmin() && p2.z() < solid_box.zmax())
    {
        return boost::optional<std::tuple<int, double, Eigen::Vector3d>>();
    }

    boost::optional<std::pair<int, double>> intersection;
    try{
        intersection  = _transform.getBlockCurrentTime(index_sequence).intersection(point,step);
    }
    catch (std::exception& e){
        std::cout << "Substrate::intersectionBlock -> " << e.what() << std::endl;
        throw std::runtime_error("Substrate::intersectionBlock -> No intersection found but segment.target() is outside block");
    }

    if (intersection)
    {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            Kernel::FT distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;
            Eigen::Vector3d normal_vector = _transform.getBlockCurrentTime(index_sequence).getNormalVector(index_face);
            double distance =  CGAL::to_double(distance_to_face);
            return std::make_tuple(index_face, distance, normal_vector);
    }
    else{
        throw std::runtime_error("WHATIS??");
    }   
}

int substrate::searchPolygon(const Eigen::Vector3d &point, int index_sequence, const std::string &frameOfReference ) const
{
    Eigen::Vector3d point_local;
    if (frameOfReference == "global")
    {
        // TODO: Check if point is out of block and if it is check that point_local is. otherwise throw error means that transform is not init properly
        transform_info transform_data = getLocalFromGlobal(point);
        point_local = transform_data.local_position;
    }
    else
    {
        point_local = point;
    }

    // TODO: Query the whole point when not using the block. Take into account when getting rid of block or implementing the transform type
    //  We will place the point in the center of the block as the tree is built with the centroid of the polygons
    Kernel::Point_3 query((double)point_local(0), (double)point_local(1), (double)point_local(2));

    Neighbor_search search(*_tree, query, 10);

    std::vector<std::pair<int, double>> indices_and_squared_distances;

    
    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        int index_polygon = _map_centroid_to_polygon.at(it->first);
        double squared_distance = CGAL::to_double(it->second); // Get the squared distance directly
        indices_and_squared_distances.push_back(std::make_pair(index_polygon, squared_distance));
    }

    
    std::sort(indices_and_squared_distances.begin(), indices_and_squared_distances.end(),
              [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              {
                  return a.second < b.second;
              });

    // Check if the point is contained in any of the sorted polygons
    for (const auto &index_squared_distance_pair : indices_and_squared_distances)
    {
        int index_polygon = index_squared_distance_pair.first;
        const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
        if (curr_poly.containsPoint(point_local))
        {
            return index_polygon;
        }
    }
    
    return -1;
}

void substrate::preComputeSubstrate(Eigen::VectorXd sequence_dt)
{
    //Compute centroid of transform block
    Eigen::Vector3d centroid = _transform.get_block_centroid();
    //Iterate over myocytes
    for (auto &myo : _myocytes)
    {
        myo.precomputePolygons(centroid, sequence_dt, _strain);
    }

    _transform.precomputeTransform(sequence_dt, _strain);
}

bool substrate::containsPoint(int index_polygon, const Eigen::Vector3d &point, int index_sequence)
{
    //return true;
    const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
    return curr_poly.containsPoint(point);
}