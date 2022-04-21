# Fitness landscape analysis on Dynamic Capacitated Arc Routing Problem. 
 
 
 ```
@inproceedings{tong2022fla,
  title={What makes the dynamic capacitated arc routing problem hard to solve: insights from fitness landscape analysis},
  author={Tong, Hao and Minku, Leandro L and Menzel, Stefan and Sendhoff, Bernhard and Yao, Xin},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference Companion},
  notes={Accepted},
  doi = {https://doi.org/10.1145/3512290.3528756}
  year={2022}
}
```

### The repository will consists of two main contents:
  - The DCARP instances used in our experiments.
  - The source code for local search and the fitness landscape analyis


### Generated Instances
- Included in the folder `instance/` 

### Sampled data for FLA
- The data size is `1.5G` so that I put it into the google drive: [link](https://drive.google.com/file/d/1zH1Xl2mavbRFSFsmHocGSuWfEkHYqeyh/view)
- Download the data from the google drive, and extract the `fla_data.zip`, and put the data into folder `analysis/`
  
### Codes
- `*.cpp` files are used for sampling and analysis
- `fla.cpp` is the main file for sampling and analysis
- You can use `Cmake` to complie the project that the `CMakeLists.txt` has been provided
- `*.py` files are used for drawing pictures

### More
Problems can be proposed in the `issus` OR you can directly contact me by the email.