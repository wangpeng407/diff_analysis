- **Description**

  - R script for differential analysis to explore markers using wilcox-test and t-test
  
  - **Example**
  
  ```
  Rscript diff_analysis.R
  ```
  
  - **Usage**
  
    Rscript diff_analysis.R mat.table group.list W prefix outdir
    
    - first argument, mat.table: column is samples, row is taxon or ko
    
    - second argument, group.list format: first column is samples, second is groups
    
    - third argument, W or T or wilcox or t.test
    
    - forth argument, prefix of outfile
    
    - fifth argument, outdir
    
  ```
  #wilcox.test()
  Rscript diff_analysis.R  input.table.xls group W prefix ./

  #t.test()
  Rscript diff_analysis.R  input.table.xls group T prefix ./
  ```

- output file:

  ID: variable ids
  
  L1_mean: mean abundance of L1 group
  
  H2_mean: mean abundance of H2 group
  
  Diff: L1_mean - H2_mean
  
  FC(fold change): L1_mean/H2_mean
  
  W-stats or T-stats: Wilcox or T-test statistics
  
  p-value and q-value
  
  
