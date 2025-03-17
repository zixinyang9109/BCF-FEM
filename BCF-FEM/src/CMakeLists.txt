add_library(bcfFEM STATIC 
       # Barycentric.cpp
        downsamplePointCloud.cpp
        
        FEM3DD.cpp
        ForceDrivenRegistration.cpp # main registration function
        Matrix3D.cpp
        pts_tools.cpp
        
        run_registration.cpp # scripts to run registrations on datasets
)