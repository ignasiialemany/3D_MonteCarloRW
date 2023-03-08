//
// Created by Alemany Juvanteny, Ignasi on 08/03/2023.
//

#include "utility_substrate.h"

std::pair<std::vector<Eigen::MatrixXd>,std::vector<Eigen::MatrixXd>> utility_substrate::read_mat_file(std::string &file)
{
    // Vertices and
    std::vector<Eigen::MatrixXd> vertices;
    std::vector<Eigen::MatrixXd> faces;    

    // Compile the full path
    std::string fullpath = "myocytes/" + file;
    const char *filepath = fullpath.c_str();

    // Read Mat file
    mat_t *geometry = Mat_Open(filepath, MAT_ACC_RDONLY);
    matvar_t *matVar, *vertex, *face;

    matVar = Mat_VarReadInfo(geometry, "myocytes");

    if (matVar){ // If the variable exist

        // Number of myocytes
        int N_m = matVar->dims[1];
        // Resizing Vertices and Faces to fit all myocytes
        vertices.resize(N_m);
        faces.resize(N_m);

        // Looping through all myocytes to retrieve vertices and faces
        for (int i = 0; i < N_m; i++){
            vertex = Mat_VarGetStructFieldByName(matVar, "Vertices", i);
            face = Mat_VarGetStructFieldByName(matVar, "Faces", i);
        
            // Reading vertices
            bool vertex_read_error = Mat_VarReadDataAll(geometry, vertex);
            if (not(vertex_read_error)){
                vertices[i] = utility_substrate::vertex_conversion(vertex);
            }
            else{
                printf("read_myocytes::unable to read vertex");
                break;
            }

            // Reading Faces
            bool face_read_error = Mat_VarReadDataAll(geometry, face);
            if (not(face_read_error)){
                faces[i] = utility_substrate::face_conversion(face);
            }
            else{
                printf("read_myocytes::unable to read face");
                break;
            }

        }
            // myocytes file has been successfully loaded
            bool flag = true;
    }
    else{
        printf("read_myocytes::non existent variable, please check your .mat file\n");
    }

    Mat_Close(geometry);
}


Eigen::MatrixXd utility_substrate::vertex_conversion(matvar_t *field){

    Eigen::MatrixXd output;
    int rows, cols;
    unsigned field_size;
    const double *Data = static_cast<const double*>(field->data) ;


    // Specify the matrix
    if (field->rank==2){
        rows = field->dims[0];
        cols = field->dims[1];
        field_size = field->nbytes / field->data_size;

        output.resize(rows, cols);
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                output(i, j) = Data[rows*j + i];
            }
        }

    }
    else{
        printf("Matrix rank not equal to 2");
    }

    return output;
}

Eigen::MatrixXd utility_substrate::face_conversion(matvar_t *field){

    Eigen::MatrixXd output;
    int rows, cols;
    unsigned field_size;
    const uint16_t *Data = static_cast<const uint16_t*>(field->data) ;


    // Specify the matrix
    if (field->rank==2){
        rows = field->dims[0];
        cols = field->dims[1];
        field_size = field->nbytes / field->data_size;

        output.resize(rows, cols);
        for (int j = 0; j < cols; j++){
            for (int i = 0; i < rows; i++){
                output(i, j) = Data[rows*j + i];
            }
        }

    }
    else{
        printf("Matrix rank not equal to 2");
    }

    return output;
}

double utility_substrate::mod(double a, double b)
{
    int r = std::fmod(a,b);
    return r < 0 ? r + b : r;
}


double utility_substrate::find_yslice(double y_global, Eigen::MatrixXd &y_slice_minmax)
{
    // To do : Use a binary search (now the worst possible method)
    bool found = false;
    for (int i = 0; i < y_slice_minmax.cols(); i++)
    {
        if (y_global > y_slice_minmax(0, i) and y_global < y_slice_minmax(1, i))
        {
            found = true;
            return y_slice_minmax(0, i);
        }
    }

    // Lost particle (Does not happen very often)
    if (not(found))
    {
        throw std::logic_error("Transform::find_yslice::where', 'Corresponding slice not found'");
        return -1;
    }
}

Eigen::Vector3d utility_substrate::rotate_y(Eigen::Vector3d &position, double theta){
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Zero(3, 3);
    rotation << std::cos(theta), 0, std::sin(theta),
                 0, 1, 0,
                -std::sin(theta), 0, std::cos(theta);
    
    return rotation*position;
}
