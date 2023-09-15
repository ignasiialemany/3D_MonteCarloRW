#ifndef OFFMESHREADER_H
#define OFFMESHREADER_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

class offMeshReader
{
public:
    offMeshReader() {}
    ~offMeshReader() {}

    bool read(const std::string& filename) {
        std::ifstream input(filename);

        if (!input.is_open()) {
            std::cerr << "Cannot open " << filename << std::endl;
            return false;
        }

        std::string header;
        input >> header;
        if (header != "OFF") {
            std::cerr << "Not a valid OFF header" << std::endl;
            return false;
        }

        int numVertices, numFaces, numEdges;
        input >> numVertices >> numFaces >> numEdges;

        _vertices.resize(numVertices, 3);
        _faces.resize(numFaces, 3); // assuming triangular mesh. Adjust if needed

        // Read vertices
        for (int i = 0; i < numVertices; ++i) {
            double x, y, z;
            input >> x >> y >> z;
            _vertices.row(i) << x, y, z;
        }

        // Read faces
        for (int i = 0; i < numFaces; ++i) {
            int n;  // Number of vertices in this face
            input >> n;
            if(n != 3) {
                std::cerr << "Non-triangular face detected. Adjust face matrix size." << std::endl;
                return false;
            }
            Eigen::Vector3d face;
            input >> face(0) >> face(1) >> face(2);
            face(0) = face(0) + 1;
            face(1) = face(1) + 1; 
            face(2) = face(2) + 1;
            _faces.row(i) = face;
        }

        return true;
    }

    const Eigen::MatrixXd getVertices() const {
        return _vertices;
    }

    const Eigen::MatrixXd getFaces() const {
        return _faces;
    }

private:
    Eigen::MatrixXd _vertices;
    Eigen::MatrixXd _faces;
};

#endif
