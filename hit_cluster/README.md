# Hit_clustering
Realiza a função de agrupar hits retornados pelo blast, se os hits possuírem coordenadas sobrepostas, impedindo que ocorra repetição de loci nas análises posteriores e aumentando o tamanho de regiões de possível homologia com o elemento utilizado como query. 

Possui a opção de realizar uma expansão e agrupamento iterativo por meio das flags
```
-ex (Size to be expanded)
-cy (Expanding cycles)
```
