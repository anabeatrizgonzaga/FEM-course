# Método dos Elementos Finitos 2D - Suavização Nodal de Tensões (Elementos T3)

**Autores:** Alexandre de Souza da Cunha e Bruna da Silva Gomes  
**Disciplina:** Método dos Elementos Finitos (MEF)  

---

## 1. Visão Geral
Este repositório contém a extensão e o aprimoramento de um código computacional escrito em Fortran para a análise de tensões e deformações em estruturas bidimensionais. O programa base utiliza o Método dos Elementos Finitos (MEF) com elementos triangulares lineares de três nós (T3), simulando cenários de **Estado Plano de Tensões (EPT)** e **Estado Plano de Deformações (EPD)**.

O objetivo principal deste trabalho foi implementar o pós-processamento de tensões para o elemento T3 (Triângulo de Deformação Constante). O código foi expandido para calcular as tensões brutas no interior dos elementos, extrapolar esses valores para os vértices geométricos e aplicar um algoritmo de suavização nodal. Isso permite a visualização de campos de tensão contínuos no software ParaView, eliminando as descontinuidades artificiais características de elementos de baixa ordem.

---

## 2. Modificações e Melhorias Implementadas

Para garantir o cálculo correto e a visualização contínua das tensões, a arquitetura original do programa foi modularizada nas seguintes frentes:

### 2.1. Isolamento da Física do Elemento (`tens3`)
A sub-rotina `tens3` foi criada do zero para processar a mecânica do elemento de forma isolada. Ela recebe os deslocamentos nodais e as matrizes cinemáticas locais, aplicando a Lei de Hooke para obter as componentes de tensão ($\sigma_{xx}$, $\sigma_{yy}$, $\sigma_{xy}$). Como o elemento T3 possui interpolação linear, a deformação resultante é constante. Logo, o algoritmo extrai um único valor de tensão posicionado no centróide do triângulo (seu único ponto de Gauss).

### 2.2. Algoritmo de Suavização Global (`pstress`)
Para resolver a descontinuidade visual gerada por tensões constantes por elemento, implementamos o motor de suavização nodal na sub-rotina `pstress`. O algoritmo opera da seguinte forma:
* **Varredura e Mapeamento:** O código mapeia a conectividade da malha, identificando quais nós pertencem a quais elementos.
* **Acúmulo Nodal:** Os valores de tensão constantes gerados pela `tens3` são acumulados diretamente nos três respectivos vértices de cada triângulo.
* **Média Simples:** O programa divide a tensão acumulada em cada nó pelo número de elementos conectados a ele, gerando um campo suavizado.
* **Cálculo de Falha:** A partir das componentes cartesianas já suavizadas, a tensão equivalente de von Mises ($\sigma_{vM}$) é computada para cada nó da malha.

### 2.3. Evolução do Pós-Processamento (`wdata` e `wvtk`)
O sistema de saída de dados foi reestruturado para suportar o armazenamento e a leitura de variáveis de tensão contínuas.
* A sub-rotina `wdata` agora imprime o bloco textual de tensões nodais suavizadas para auditoria e verificação analítica.
* A sub-rotina `wvtk` teve sua sintaxe corrigida e expandida. As tensões foram associadas aos nós geométricos da malha utilizando a tag `POINT_DATA` e exportadas como `SCALARS`, garantindo que o ParaView renderize um gradiente de cores contínuo ao longo da estrutura.

---

## 3. Metodologia de Validação e Testes

Para atestar a precisão da formulação matemática e o funcionamento da rotina de suavização, estruturamos os seguintes ensaios numéricos:

### Ensaio 1: Tração Pura em Viga Engastada (Validação de Código)
* **Objetivo:** Provar a exatidão das rotinas `tens3` e `pstress` em um cenário onde não há variações de momento fletor.
* **A Simulação:** Uma viga foi engastada em uma das extremidades e submetida a uma carga de tração uniforme na extremidade oposta.
* **Resultado Obtido:** O elemento T3 reproduziu o estado de tração pura exatamente, demonstrando que o campo de tensão constante pertence ao espaço de aproximação do elemento. O programa reportou um desvio de precisão praticamente nulo (ruído numérico na casa de $10^{-5}$ para as tensões transversais), validando toda a álgebra matricial do código.

### Ensaio 2: Flexão Simples (Análise de Limitações Inerentes)
* **Objetivo:** Avaliar o comportamento numérico do elemento e do algoritmo de suavização sob a presença de gradientes acentuados de tensão.
* **A Simulação:** Uma viga biapoiada submetida a um carregamento distribuído constante ao longo de seu comprimento. 
* **Resultado Obtido:** A simulação obteve uma tensão normal máxima abaixo da solução analítica exata. O teste confirmou a limitação intrínseca da formulação do T3: a rigidez matemática excessiva (cisalhamento espúrio) penaliza a deformação em curvatura. Isso comprovou que, embora o código esteja correto e a suavização seja eficiente para a continuidade visual, a baixa ordem do elemento limita a captura precisa da física da flexão.

---

## 4. Estrutura do Código-Fonte

O repositório contém a versão consolidada em Fortran:

* `mef.f`: Arquivo principal do programa contendo a rotina controladora (`contr`), solver PCG, formulação das matrizes de rigidez (`elmt01`, `elmt02`) e as novas sub-rotinas responsáveis pela física e pós-processamento das tensões (`pstress`, `tens3`, `wdata`, `wvtk`).
