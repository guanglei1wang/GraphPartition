# Users' Manual

Authors: Guanglei Wang and Hassan Hijazi 
Affiliations: The Australian National University, ACTON 2601, Canberra, Australia; Los Alamos National Laboratory, New Mexico, USA.

## Introduction
The executable `min_k_part` has been compiled and tested on a MPC server. It is under the folder `/home/102017-00081/MKP/bin/`.  

Users can use it to reproduce numerical results reported in the associated manuscript "Exploiting sparsity for the min k-partition problem" (\#102017-00081). 

The Syntax section explains how the code works; the Data section provides information on the content and the location of each set of the test instances.

## Syntax
The executable `min_k_part` solves min k-partition problem.
It takes four ordered input arguments, namely:

1. the input graph data `*.txt`
2. the integer parameter on the number of partitions `k` (k > 2, we tested k = 3, 4 in our numerical tests) 
3. the formulation name (`MIP` or `Node_edge` or `MIP_tree` or `SDP` or `SDP_tree`) 
4. the boolean type for continuous relaxation (`true` for relax, `false` for integer program). 

Note that `MIP` corresponds to Model 1 in the paper, `Node_edge` to Model 2, `SDP` to Model 3, `MIP_tree` to Model 4, `SDP_tree` to Model 5. 


## Data 

All the numerical tests results will also be outputted to the file `MkP_result.txt`. 

* All test instances mentioned in the paper can be found in the folder `~/MKP/data_sets/Minkcut/` and they are in `.txt` format.  

* For each problem instance, they are generally named by `graphtype + '_' +the number of vertices in the graph`. For instance, `spinglass2g_1010.txt` refers to the `spinglass2g` problem instance with `10 x 10` vertices.

	* As stated at page 14 of the manuscript. The sparsity of the `band` graphs is parametrized by `k`. So we additionally add `_3` or `_4` to specify this parameter. For instance, `band100_3` refers to the band graph with 100 nodes and edge set {(i ,j) : i - j <= 4} (i.e, 390 edges). 
	* For random graphs, `random10_150` refers to a random graph with 150 nodes and edge density around 10%. 


## Illustration

#### Example 1: ILP formulations
Assume that you are in the folder `MKP/bin/` and want to valid the result of `spinglass2g` with `k = 3` and `|V|= 6 x 6`in Table 1, you need to type the following to get results for Model 1

`./min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP true` 

`./min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP false`

 
 and  
  
 `./min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP_tree true`

 `./min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP_tree false`

to get results for Model 4. 

You should get similar results as below in `MKP/build/MkP_result.txt`:
> 3, 36, 12.1644,  -2.14607e+06  
> 3, 36, 12.9085,  -2.14607e+06  
> 3, 36, 0.159482, -2.14607e+06  
> 3, 36, 0.368595, -2.14607e+06

Each row shows the corresponding result of each command line above. Thus we have four rows. The first column is parameter k; the second is the number of nodes which is `6x6=36`; the third column is the CPU time in seconds; the fourth is the optimal value, by which you can calculate the optimality gap (0 for this test instance) for their continuous relaxations. 

Since the data in `MkP_result.txt`is separated by comma,  you can open this `csv (comma-separated-value)` file using MS Excel and calculate the optimality gap easily. 

#### Example 2: ISDP formulations
We illustrate tests of ISDP formulations (Model 3 and Model 5 in our paper) with instance `(4, 11x11), spinglass2g` in Table 3. You need to type the following to get results for Model 3:

`./min_k_part ../data_sets/Minkcut/spinglass2g_1111.txt 4 SDP true`


and 
`./min_k_part ../data_sets/Minkcut/spinglass2g_1111.txt 4 SDP_tree true`

for Model 5. 

You should get similar results as below in `MKP/build/MkP_result.txt`:

> 4,121,145.845, -8.44258e+06  
> 4,121,1.12703, -8.45322e+06

### Run multiple tests using the bash script
If you want to run multiple the problem instances in the paper, you can use the bash script called `exp.sh` that we provided in folder `MKP`. 

1. Copy the `exp.sh` file to your build folder. 

2. Modify the content of the `exp.sh` file to run corresponding set of numerical experiments. 

3. In your build folder, type `./exp.sh`. 

4. check the result in `MkP_result.txt`. 


## Contact 

Guanglei Wang: <guanglei.wang@anu.edu.au>

Hassan Hijazi: <hassan.hijazi@anu.edu.au>


