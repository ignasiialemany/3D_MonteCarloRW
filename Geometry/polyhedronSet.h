#ifndef POLYHEDRONSET_H
#define POLYHEDRONSET_H

#include "../Geometry/polygon.h"
#include "../Utility/offMeshReader.h"

class polyhedronSet
{
public:
    polyhedronSet(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input);
    polyhedronSet(const std::string &filename);
    polyhedronSet()=default;
    ~polyhedronSet()=default;

    const polygon& getPolygon(int index)const {return _polygonSet[index];};
    void precomputePolygons(Eigen::Vector3d &centroid, Eigen::VectorXd &sequence_dt, std::function<double(double)> strain = [](double){return 0.;});
    Eigen::Vector3d getUnstrainedCentroid(){return _unstrainedcentroid;};
    std::vector<polygon> _polygonSet;


private:

    int _ntimeSteps;
    std::function<double(double)> _strain;
    Eigen::MatrixXd _unstrainedvertices;
    Eigen::MatrixXd _faces;
    Eigen::Vector3d _unstrainedcentroid;

    Eigen::MatrixXd compute_vertex(double strain, Eigen::Vector3d &centroid_strain);






};

#endif