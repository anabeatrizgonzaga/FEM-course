# Implementação de Cargas de Corpo (Peso Próprio) em Elementos T3 e Q4

# Descrição

Este trabalho consiste na implementação de cargas de corpo (peso próprio) em um código de Método dos Elementos Finitos (MEF) para análise bidimensional de elasticidade.

A implementação permite que forças distribuídas de origem gravitacional sejam convertidas em forças nodais equivalentes e incorporadas automaticamente ao vetor global de carregamentos durante a montagem do sistema de equações.

Foi desenvolvida a sub-rotina `bodyforce`, responsável por calcular as forças nodais equivalentes a partir da densidade do material e das acelerações aplicadas nas direções x e y. A rotina foi integrada ao processo de montagem dos elementos, permitindo que tanto elementos triangulares de três nós (T3) quanto elementos quadrilaterais de quatro nós (Q4) considerem carregamentos de corpo de forma transparente ao usuário.

O código também realiza a exportação dos resultados para o formato VTK, permitindo a visualização da malha, dos deslocamentos e da configuração deformada em softwares de pós-processamento, como o ParaView.


# Formulação Teórica

As cargas de corpo representam forças distribuídas atuando continuamente sobre o volume do material, como o peso próprio decorrente da ação da gravidade.

A carga de corpo por unidade de volume é definida por

\[
\mathbf{b}=\rho\mathbf{g}
\]

onde:

- **ρ** é a densidade do material;
- **g** = (gx, gy) representa as acelerações aplicadas nas direções x e y.

O vetor de forças nodais equivalentes é obtido por

\[
\mathbf{f}_b=\int_{\Omega_e}\mathbf{N}^T\mathbf{b}\,t\,d\Omega
\]

em que:

- **N** é a matriz das funções de forma;
- **t** é a espessura do elemento;
- **Ωe** representa o domínio do elemento.

## Elemento T3

Para elementos triangulares lineares (T3), as funções de forma são lineares e a carga de corpo é constante. Dessa forma, a integração é realizada analiticamente, resultando em

\[
\mathbf{f}_b=
\frac{\rho At}{3}
\begin{Bmatrix}
g_x\\
g_y\\
g_x\\
g_y\\
g_x\\
g_y
\end{Bmatrix}
\]

onde:

- **A** é a área do elemento.

Assim, o peso próprio do elemento é distribuído igualmente entre os seus três nós.

## Elemento Q4

Para elementos quadrilaterais bilineares (Q4), o vetor de forças nodais equivalentes é obtido por integração numérica utilizando quadratura de Gauss 2×2:

\[
\mathbf{f}_b=
\sum_{i=1}^{4}
\mathbf{N}^T(\xi_i,\eta_i)
\rho\mathbf{g}
t
|J(\xi_i,\eta_i)|
w_i
\]

em que:

- **(ξi, ηi)** são os pontos de integração;
- **J** é o Jacobiano da transformação isoparamétrica;
- **wi** são os pesos da quadratura de Gauss.


# Como Compilar e Executar

1. Abra a solução no Microsoft Visual Studio.
2. Compile o projeto.
3. Execute o programa.
4. Informe:
   - arquivo de entrada (`.dat`);
   - arquivo de saída (`.out`);
   - arquivo VTK (`.vtk`).
5. Após a execução, visualize o arquivo `.vtk` no ParaView.


# Exemplo

Foram utilizados três casos de teste para verificar a implementação das cargas de corpo em elementos T3 e Q4, avaliando tanto problemas simples quanto geometrias mais complexas.

## Exemplo 1 – Placa discretizada com elementos T3 (`Ex1_tria_200e.dat`)

Modelo bidimensional de uma placa discretizada por elementos triangulares lineares (T3).

**Objetivo**

Verificar o cálculo das forças nodais equivalentes devido às cargas de corpo utilizando elementos triangulares, comparando os deslocamentos e a distribuição das forças com os resultados esperados.

---

## Exemplo 2 – Placa discretizada com elementos Q4 (`Ex1_quad_200e.dat`)

Modelo equivalente ao Exemplo 1, porém discretizado com elementos quadrilaterais bilineares (Q4).

**Objetivo**

Validar a implementação para elementos quadrilaterais, verificando a integração numérica por quadratura de Gauss 2×2 e comparando os resultados com aqueles obtidos para a malha triangular.

---

## Exemplo 3 – Corpo de prova tipo Dogbone (`dogbone.dat`)

Modelo bidimensional de um corpo de prova do tipo dogbone, discretizado por elementos quadrilaterais.

**Objetivo**

Avaliar o desempenho da implementação em uma geometria mais complexa, verificando a correta montagem das cargas de corpo, a obtenção dos deslocamentos nodais e a geração dos arquivos VTK para visualização dos resultados no ParaView.

---

## Resultados obtidos

Os três casos de teste apresentaram resultados consistentes, demonstrando o correto funcionamento da implementação. Em todos os exemplos observou-se:

- cálculo correto das forças nodais equivalentes devido às cargas de corpo;
- montagem consistente do vetor global de carregamentos;
- convergência do solucionador;
- obtenção dos deslocamentos nodais esperados;
- geração automática dos arquivos VTK para pós-processamento;
- visualização da malha deformada e do campo de deslocamentos no ParaView.

# Verificação

A implementação foi validada por meio dos seguintes procedimentos:

# Verificação

A implementação foi validada por meio dos seguintes procedimentos:

- comparação da solução analítica do elemento T3 com as forças nodais equivalentes obtidas numericamente;
- verificação da conservação da força total aplicada, comparando o peso próprio do elemento com a soma das forças nodais equivalentes;
- testes utilizando acelerações independentes nas direções x e y;
- comparação dos resultados obtidos para elementos T3 e Q4 em problemas equivalentes;
- inspeção visual dos deslocamentos e da configuração deformada utilizando o ParaView;
- verificação da convergência do solucionador para todos os casos analisados.


# Referências

1. RIBEIRO, F. L. B. **Introdução ao Método dos Elementos Finitos**. COPPE/UFRJ, Rio de Janeiro, 2004.

2. ZIENKIEWICZ, O. C.; TAYLOR, R. L.; ZHU, J. Z. **The Finite Element Method: Its Basis and Fundamentals**. 7th ed. Butterworth-Heinemann, 2013.

3. LOGAN, D. L. **A First Course in the Finite Element Method**. 6th ed. Cengage Learning, 2017.