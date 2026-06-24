# Hit_clustering
Realiza a função de agrupar hits retornados pelo blast, se os hits possuírem coordenadas sobrepostas, impedindo que ocorra repetição de loci nas análises posteriores e aumentando o tamanho de regiões de possível homologia com o elemento utilizado como query. 

Possui a opção de realizar expansão e agrupamento iterativo por meio das flags, dessa forma é possível é possível agrupar, não apenas por sobreposição,\n
mas também por proximidade em pares de base entre os hits (flag -ex).

```
-ex (Tamanho em bp em que hits próximos serão agrupados)
-cy (ciclos de expansão)
```
Além de realizar a busca via blast e expansão, o programa retorna a sequência em FASTA das regiões expandidas utilizando o comando bedtools get-fasta.\n 
Além disso sequências FASTA extraídas são agrupadas utilizando o programa CD-HIT. 

Demais flags:

```
-o (Output file name)
-s (Separador utilizado nos arquivo de hits, padrão do programa blastn é \t)

```
