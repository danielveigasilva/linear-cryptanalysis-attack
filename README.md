# linear-cryptanalysis-attack

Esta aplicação realiza um ataque de criptoanálise linear a partir de uma S-BOX 256, por meio de textos planos e cifrados conhecidos calculamos a chave privada.


A aplicação gera dois arquivos de saída: 

* biastable.csv: Contendo a tabela de bias gerada após a análise da S-BOX
* relatorio.tex: Relatório em LATEX contendo a tabela de bias (resumida caso seja muito grande), as equações utilizadas, suas probabilidades, quantidade de interações e a chave encontrada.

Para a execução utilize o seguinte comando:

```console
gcc linear.c && ./a.out && pdflatex relatorio.tex
```

Caso não possua o LATEX instalado:

```console
gcc linear.c && ./a.out
```
