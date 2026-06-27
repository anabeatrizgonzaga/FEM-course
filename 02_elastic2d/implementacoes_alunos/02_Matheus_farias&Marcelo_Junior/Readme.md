# Método dos Elementos Finitos 2D - Análise Estrutural (Elementos Q4)

**Autores:** Matheus Sousa Farias e Marcelo Junior  
**Disciplina:** Método dos Elementos Finitos (MEF)  

---

## 1. Visão Geral do Projeto

Este repositório contém a extensão e o aprimoramento de um código computacional escrito em Fortran para a análise de tensões e deformações em estruturas bidimensionais. O programa utiliza o Método dos Elementos Finitos (MEF) com elementos triangulares de três nós (T3) e quadrilaterais de 4 nós (Q4), sendo capaz de simular cenários de **Estado Plano de Tensões (EPT)** e **Estado Plano de Deformações (EPD)**.

O objetivo principal das implementações deste trabalho foi criar uma subrotina para extrair os resultados de deslocamento e transformá-los em campos de tensões, permitindo a visualização gráfica das estruturas submetidas a carregamentos. 
Isto é, criar uma subrotina que recupera σx, σy, τxy a partir dos deslocamentos nodais para os elementos Q4 (elmt03 e elmt04), recuperar σx, σy, τxy nos pontos de Gauss 2×2 e imprimir as saídas de tensões nos pontos de graus no arquivo de saída, em formato .vtk

---

## 2. Modificações e Melhorias Implementadas

Para que o cálculo das tensões fosse realizado, o código original foi expandido nas seguintes partes:

### 2.1. Diferenciação da Física (EPT e EPD)
A Matriz Constitutiva, que dita o comportamento elástico do material, foi ajustada para identificar automaticamente qual formulação física utilizar com base nos dados de entrada. Isso garante que o efeito de compressão ou expansão lateral (Efeito de Poisson) seja calculado corretamente, respeitando a espessura da estrutura física.

### 2.2. Integração nos Pontos de Gauss
Na mecânica dos Elementos Finitos, as tensões não são calculadas diretamente nos vértices (nós), mas sim em pontos de integração numéricos no interior de cada elemento geométrico, conhecidos como Pontos de Gauss. Implementamos a rotina matemática para avaliar as matrizes cinemáticas nestes pontos internos, garantindo a obtenção das tensões (σx, σy, τxy) no núcleo de cada quadrado.

### 2.3. Evolução do Pós-Processamento
Para que o software de visualização (ParaView) consiga plotar o mapa de cores das tensões, precisamos exportar os resultados matemáticos de forma legível. Isso foi feito em duas etapas evolutivas, gerando duas versões distintas do código:

**A. Primeira Etapa: A Abordagem por Média (Discreta)**
Na primeira versão (`Solucao-mef_media`), o código pega os valores de tensão calculados nos pontos internos de Gauss e extrai uma média simples para representar o elemento inteiro. Como resultado, cada quadrado na malha recebe um valor único e uma cor sólida. O mapa visual gerado é discreto (semelhante a blocos ou pixels). É útil para verificação global do equilíbrio, mas possui a limitação de não refletir a transição contínua de forças que ocorre.

**B. Segunda Etapa: A Abordagem Suavizada (Projeção Nodal)**
Para resolver a descontinuidade visual e simular o comportamento contínuo dos materiais reais, desenvolvemos a versão suavizada (`Solucao-mef_suavizado`). O processo ocorre em dois passos:
1. **Extrapolação:** O código pega as tensões exatas no interior do elemento e as projeta matematicamente de volta para os seus quatro vértices (nós).
2. **Média Nodal:** Como os elementos formam uma malha conectada, um mesmo nó pode pertencer a até quatro quadrados diferentes, recebendo valores de tensão ligeiramente distintos de cada vizinho. O programa agrupa todas essas contribuições que chegam a um mesmo nó e calcula uma média local. 

O resultado final é um único valor de tensão por nó. Isso permite que o software de pós-processamento crie interpolações entre os nós, gerando um mapa em gradiente contínuo, sem quebras bruscas de cor.

---

## 3. Metodologia de Validação e Testes

Para provar a precisão das matrizes implementadas, criamos cinco cenários de testes, disponíveis na pasta de exemplos:

### Exemplo 1: Malha de 25 Elementos (Teste de Sintaxe e Conectividade)
* **Objetivo:** Calibrar o código para garantir a correta leitura, alocação de memória e escrita do formato estruturado do arquivo VTK.
* **Configuração:** Uma placa bidimensional subdividida em 25 elementos Q4, submetida a um carregamento genérico. 
* **Resultado:** O teste confirmou que as rotinas de numeração de nós e conectividade dos elementos estavam perfeitamente sincronizadas, gerando a geometria correta no software de pós-processamento.

### Exemplo 2: Tração Pura em 1 Elemento (Usado para validação)
* **Objetivo:** Provar a exatidão da Matriz Jacobiana e do motor de cálculo de tensões do código Fortran, isolando a física do problema para um cenário mais simples e previsível.
* **A Simulação:** Modelamos um único quadrado perfeito com dimensões de 1 x 1 metro. Atribuímos a ele um material teórico com Módulo de Elasticidade (E) de 1000 e Coeficiente de Poisson zerado (para impedir interferências tridimensionais).
* **Condições de Contorno e Carga:** Travamos completamente os dois nós da face esquerda (impedindo qualquer movimentação em X e Y, simulando um engaste em uma parede). Na face direita, aplicamos uma força tracionadora total de 100 unidades (dividida em 50 para cada nó), puxando o quadrado para fora.
* **O que era esperado (Cálculo Analítico):**
  * Pela física básica, a Tensão no eixo X deve ser igual à Força dividida pela Área (Tensão = F/A). Sendo F = 100 e A = 1, a Tensão X esperada era exatamente 100.
  * Como a força era perfeitamente horizontal, as Tensões em Y e as Tensões de Cisalhamento (XY) deveriam ser exatamente zero.
  * Pela Lei de Hooke (Tensão = Módulo de Elasticidade * Deformação), o deslocamento em X do lado direito deveria ser exatamente 0.1.
* **Resultado Obtido:** O arquivo gerado pelo nosso código acusou, em todo o elemento, o valor de `1.0e+02` (100) para a Tensão X, exatos `0.1` para o deslocamento, e valores nulos para os demais esforços. Isso comprova que toda a matemática das rotinas implementadas está correta e funcional.

### Exemplo 3: Malha de 100 Elementos com Efeito de Poisson (Validação da Suavização)
* **Objetivo:** Testar o código sob um estresse matemático maior e validar o algoritmo de suavização nodal (projeção contínua de cores) no ParaView.
* **A Simulação:** Uma placa refinada com 100 elementos e 121 nós, engastada em toda a sua lateral esquerda. Desta vez, introduzimos um Coeficiente de Poisson real (0.3), forçando o programa a acoplar o cálculo dos eixos X e Y.
* **Carga e Resultado:** Aplicamos um carregamento vertical de cima para baixo no topo da estrutura. Como resultado, o mapa do ParaView gerou campos de tensão contínuos e sem quebras abruptas de cor, exibindo uma clara zona de tração na parte superior próxima ao engaste e uma zona de compressão na parte inferior, retratando o comportamento real de flexão da estrutura.

### Exemplo 4: Compressão Confinada em 1 Elemento (Validação EPD)
* **Objetivo:** Comprovar a resposta do código à formulação de Estado Plano de Deformações (EPD - `elmt03`), validando o cálculo transversal gerado pelo Efeito de Poisson sob confinamento total.
* **A Simulação:** Um bloco de 1 x 1 metro (E=1000, Poisson=0.3), totalmente travado na base e com as laterais impedidas de expandir no eixo X. Uma carga de compressão vertical de -100 foi aplicada no topo.
* **O que era esperado (Cálculo Analítico):**
  * A Tensão Y esperada era exatamente -100.
  * Como o elemento EPD simula uma estrutura infinitamente profunda (travada no eixo Z), o estufamento lateral fica duplamente confinado. Matematicamente, a Tensão X gerada pela parede deveria ser exatamente -42.85.
  * O deslocamento vertical no topo deveria ser -0.074, refletindo a alta rigidez da estrutura confinada.
* **Resultado Obtido:** Os resultados do ParaView e do arquivo de saída cravaram exatamente nos valores analíticos (-100 em Y, -42.85 em X e -0.074 de deslocamento), comprovando a precisão da Matriz Constitutiva para EPD.

### Exemplo 5: Compressão Confinada em 1 Elemento (Validação EPT)
* **Objetivo:** Provar a diferenciação física do código ao trocar para o Estado Plano de Tensões (EPT - `elmt04`), simulando o mesmo carregamento em uma chapa fina.
* **A Simulação:** O exato mesmo modelo do Exemplo 4 (bloco 1 x 1, carga de -100, travamento lateral em X), mas utilizando o elemento tipo 4 e declarando a espessura da chapa como 1.0.
* **O que era esperado (Cálculo Analítico):**
  * A Tensão Y continuaria sendo -100.
  * Como uma chapa fina (EPT) está livre para "estufar" na direção Z, a pressão exercida contra as paredes laterais (eixo X) é menor. A Tensão X esperada cairia para exatamente -30.
  * Por estar livre em Z, a estrutura se torna menos rígida que no EPD. O deslocamento vertical esperado era maior: -0.091.
* **Resultado Obtido:** O programa registrou perfeitamente a queda de rigidez, retornando a Tensão X de -30 e deslocamento de -0.091. Este teste cruzado serve como prova de que o software interpreta corretamente a mecânica tridimensional por trás de cada formulação.

---

## 4. Estrutura de Diretórios

O repositório foi organizado para separar de forma limpa o código-fonte dos dados de validação gráfica:

* `scr/`: Contém os códigos em Fortran.
    * `Solucao-mef_media.f` (Gera as tensões médias por elemento).
    * `Solucao-mef_suavizado.f` (Gera o campo de tensões nodais contínuas).
* `exemplos/`: Diretório de testes.
    * `Media/`: Contém os cinco exemplos calculados com a primeira abordagem.
    * `Suavizado/`: Contém os cinco exemplos calculados com a projeção nodal avançada.
    * *Nota:* Dentro de cada diretório de exemplo estão os arquivos de entrada de dados (`.dat`), os relatórios brutos de saída, o arquivo interpretável (`.vtk`) e uma pasta dedicada às imagens de comprovação do ParaView.