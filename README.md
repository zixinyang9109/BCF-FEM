# Boundary Constraint-Free Biomechanical Model-Based Surface Matching for Intraoperative Liver Deformation Correction

## Datasets

Datasets include our constructed phantoms, PBSM dataset, and [Sparse data challenge dataset](https://osf.io/etpq9).

The preprocessed dataset and our raw results can be downloaded via:

- [Demo](https://drive.google.com/drive/folders/1eiGozqjyq-7D0obMpAE-W-mebk_UKNp8?usp=sharing) 
- [Our phatnoms](https://drive.google.com/drive/folders/1rrtw6r6oDdFevubUwgiszS884Aef9xrc?usp=sharing)
- [PBSM dataest](https://drive.google.com/drive/folders/1QAGupR4feKBT0FdsBJXGhF74WL53JCmk?usp=sharing) 
- [Sparse dataest](https://drive.google.com/drive/folders/1TzAVQzJkgcccABnnGOT5Wfy1Rp8xA0Y3?usp=sharing)

## Our Code

### Run Non-rigid Organ Registration
Our registration code is within the BCF-FEM folder. 

The only dependence is the [eigen library](https://drive.google.com/drive/folders/1TnSgEn-Km1tiKHBdASub9rFEpTOj6R8J?usp=sharing). Download and put it inside the BCF-FEM folder, or change its path in the CMakeLists.txt. 

The code is developed on Windows 10 using [Visual Studio 2019](https://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=msvc-170). You can run the following files within the VS or use built .exe files. Please set your EXECUTABLE_OUTPUT_PATH and LIBRARY_OUTPUT_PATH in the CMakeLists.txt. 



```
# Please change the data paths to your paths 
app/Demo.cpp                       
app/run_our_phantom.cpp
app/run_pbsm_dataset.cpp
app/run_sparse_dataset.cpp
```

### Evaluation

```
# change the data paths to yours
Evaluation/eval_our_dataset.py
Evaluation/eval_pbsm_dataset.py
Evaluation/eval_sparse_dataset.py
```

### Test on Other Datasets

We use the [.off](https://en.wikipedia.org/wiki/OFF_(file_format)) format to store data. If you want to transfer other formats such as .vtk and .ply to .off. Please see the functions within the
```
Data_Format_Prepare/format_transfer.py
```
## Reference
The reference for our phantoms:
```bibtex
@article{yang2024boundary,
  title={Boundary Constraint-free Biomechanical Model-Based Surface Matching for Intraoperative Liver Deformation Correction},
  author={Yang, Zixin and Simon, Richard and Merrell, Kelly and Linte, Cristian A},
  journal={IEEE Transactions on Medical Imaging},
  year={2024},
  publisher={IEEE}
}
```

Three references for the Sparse Data Challenge dataset:
```bibtex
@inproceedings{brewer2019image,
  title={The image-to-physical liver registration sparse data challenge},
  author={Brewer, E Lee and Clements, Logan W and Collins, Jarrod A and Doss, Derek J and Heiselman, Jon S and Miga, Michael I and Pavas, Chris D and Wisdom III, Edward H},
  booktitle={Medical Imaging 2019: Image-Guided Procedures, Robotic Interventions, and Modeling},
  volume={10951},
  pages={364--370},
  year={2019},
  organization={SPIE}
}

@article{collins2017improving,
  title={Improving registration robustness for image-guided liver surgery in a novel human-to-phantom data framework},
  author={Collins, Jarrod A and Weis, Jared A and Heiselman, Jon S and Clements, Logan W and Simpson, Amber L and Jarnagin, William R and Miga, Michael I},
  journal={IEEE transactions on medical imaging},
  volume={36},
  number={7},
  pages={1502--1510},
  year={2017},
  publisher={IEEE}
}

@article{heiselman2024image,
  title={The image-to-physical liver registration sparse data challenge: comparison of state-of-the-art using a common dataset},
  author={Heiselman, Jon S and Collins, Jarrod A and Ringel, Morgan J and Peter Kingham, T and Jarnagin, William R and Miga, Michael I},
  journal={Journal of Medical Imaging},
  volume={11},
  number={1},
  pages={015001--015001},
  year={2024},
  publisher={Society of Photo-Optical Instrumentation Engineers}
}
```

The reference for PBSM dataset:
```bibtex
@article{suwelack2014physics,
  title={Physics-based shape matching for intraoperative image guidance},
  author={Suwelack, Stefan and R{\"o}hl, Sebastian and Bodenstedt, Sebastian and Reichard, Daniel and Dillmann, R{\"u}diger and dos Santos, Thiago and Maier-Hein, Lena and Wagner, Martin and W{\"u}nscher, Josephine and Kenngott, Hannes and others},
  journal={Medical physics},
  volume={41},
  number={11},
  pages={111901},
  year={2014},
  publisher={Wiley Online Library}
}
```


## Related Software/ Projects

- [MeshLab](https://www.meshlab.net/) (Visualization) 
- [PyVista](https://pyvista.org/) (Data processing, evaluation, and visualization)
- [Gmsh](https://gmsh.info/) (Volumetric mesh generation)
- [GMM-FEM](https://github.com/siavashk/GMM-FEM) (Baseline comparison)
- [V2S-Net](https://gitlab.com/nct_tso_public/Volume2SurfaceCNN) (Baseline comparison)
- [Adjoint-elastic-registration](https://github.com/gmestdagh/adjoint-elastic-registration) (Data term formulation)

