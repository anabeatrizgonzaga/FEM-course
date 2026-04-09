# MEF — Programa de Elementos Finitos para Elasticidade Plana

> Código educacional em FORTRAN 77 para análise estrutural 2D por elementos finitos.  
> Disciplina de Elementos Finitos — Engenharia Civil

---

## Sumário

1. [O que este programa faz](#1-o-que-este-programa-faz)
2. [Como compilar e executar](#2-como-compilar-e-executar)
3. [Fundamentos teóricos](#3-fundamentos-teóricos)
4. [Arquitetura do programa](#4-arquitetura-do-programa)
5. [Tipos de elementos disponíveis](#5-tipos-de-elementos-disponíveis)
6. [Formato do arquivo de entrada `.dat`](#6-formato-do-arquivo-de-entrada-dat)
7. [Exemplos resolvidos](#7-exemplos-resolvidos)
8. [Formato do arquivo de saída](#8-formato-do-arquivo-de-saída)
9. [O solver PCG explicado](#9-o-solver-pcg-explicado)
10. [Armazenamento Skyline explicado](#10-armazenamento-skyline-explicado)
11. [Referências](#11-referências)

---

## 1. O que este programa faz

Este programa resolve problemas de **elasticidade linear plana** pelo Método dos Elementos Finitos (MEF). Dado um sólido bidimensional com geometria, material, apoios e cargas definidos, ele calcula os **deslocamentos nodais** em cada ponto da malha.

O problema central que o programa resolve é o sistema linear:

```
K · u = f
```

onde:
- **K** é a matriz de rigidez global (propriedades do material + geometria da malha)
- **u** é o vetor de deslocamentos nodais (incógnitas)
- **f** é o vetor de forças nodais equivalentes (cargas aplicadas)

O programa suporta dois casos físicos e dois tipos de elemento, combinados em quatro opções:

| Código | Elemento | Hipótese física | Parâmetros |
|:------:|----------|-----------------|------------|
| `1` | Triângulo T3 | Estado Plano de **Deformações** (EPD) | E, ν |
| `2` | Triângulo T3 | Estado Plano de **Tensões** (EPT) | E, ν, espessura |
| `3` | Quadrilátero Q4 | Estado Plano de **Deformações** (EPD) | E, ν |
| `4` | Quadrilátero Q4 | Estado Plano de **Tensões** (EPT) | E, ν, espessura |

> **EPD (Estado Plano de Deformações):** deformação fora do plano é nula (ε_z = 0). Aplicável a barragens, seções de túneis, fundações — estruturas longas na direção z.
>
> **EPT (Estado Plano de Tensões):** tensão fora do plano é nula (σ_z = 0). Aplicável a chapas finas carregadas no plano.

---

## 2. Como compilar e executar

### Compilação

```bash
# Com gfortran:
gfortran -O2 -o mef program.for

# Com ifort (Intel):
ifort -O2 -o mef program.for
```

### Execução

```bash
./mef
```

O programa solicita dois nomes de arquivo interativamente:

```
Arquivo de dados:   Ex1.dat       <- arquivo de entrada que você prepara
Arquivo de saida:   resultado.out <- arquivo de resultados gerado pelo programa
```

### Saída esperada no terminal

```
Numero de equacoes (graus de liberdade livres): 220
Tamanho do perfil skyline (ncs): 3540
Iniciando solver PCG...
PCG convergiu em 87 iteracoes.
Numero de equacoes    :    220
Norma de energia      :   0.12345678
Calculo concluido. Resultados gravados com sucesso.
```

---

## 3. Fundamentos teóricos

### 3.1 A ideia central do MEF

O MEF divide o domínio (a estrutura) em pequenas regiões chamadas **elementos**. Dentro de cada elemento, o campo de deslocamentos é aproximado por **funções de interpolação** (também chamadas funções de forma), que são polinômios simples.

```
Domínio real              Malha de elementos
  ___________               ___________
 /           \             /\  /\  /\ /\
/    sólido   \    -->    /  \/  \/  X  \
\             /           \  /\  /\ / \  /
 \___________/             \/  \/  /   \/
```

A solução em todo o domínio é então determinada calculando os deslocamentos nos **nós** (pontos de intersecção dos elementos).

### 3.2 Da física à matriz de rigidez

O ponto de partida é o **Princípio dos Trabalhos Virtuais**. Para um sólido em equilíbrio, qualquer deslocamento virtual δu satisfaz:

```
INT INT (δε : σ) dΩ  =  INT INT (δu · b) dΩ  +  INT (δu · t) dΓ
```

Substituindo a relação constitutiva **σ = D ε** e a aproximação **u ≈ N û** (onde N são as funções de forma e û os deslocamentos nodais), chega-se ao sistema:

```
K · u = f
```

A **matriz de rigidez elementar** é calculada como:

```
Ke = INT INT  B^T · D · B  dΩ
```

onde:
- **B** = matriz de deformação (derivadas das funções de forma N em relação a x e y)
- **D** = matriz constitutiva (propriedades do material)

### 3.3 Matrizes constitutivas D

**Estado Plano de Deformações (EPD):**

```
         E(1-v)         |  1        v/(1-v)      0           |
D = ─────────────────  ·|  v/(1-v)     1         0           |
     (1+v)(1-2v)        |  0            0    (1-2v)/(2(1-v)) |
```

**Estado Plano de Tensões (EPT):**

```
      E        |  1    v      0      |
D = ──────── · |  v    1      0      |
    (1 - v²)   |  0    0   (1-v)/2   |
```

### 3.4 Integração numérica de Gauss para o Q4

Para o quadrilátero, a integral de Ke não tem solução analítica simples, então usa-se **Quadratura de Gauss 2x2** (4 pontos):

```
INT INT f(r,s) dr ds  ≈  SUM wi · wj · f(ri, sj)     (wi = wj = 1)
```

Os 4 pontos de Gauss têm coordenadas (r, s) = (+/- 1/sqrt(3), +/- 1/sqrt(3)).

```
     s
     |
  3--+--4      Coordenadas:
  |  |  |      Ponto 1: (+1/sqrt(3), -1/sqrt(3))
--+--+--+-- r  Ponto 2: (+1/sqrt(3), +1/sqrt(3))
  |  |  |      Ponto 3: (-1/sqrt(3), +1/sqrt(3))
  1--+--2      Ponto 4: (-1/sqrt(3), -1/sqrt(3))
     |
```

O valor `1/sqrt(3) = 0.577350269189626` aparece explicitamente no código.

---

## 4. Arquitetura do programa

O programa é organizado em sub-rotinas com responsabilidades bem definidas. O fluxo completo é:

```
PROGRAMA PRINCIPAL (mef)
|
+-- contr()   controlador geral
    |
    +-- [1] rdata()    Le arquivo .dat (nos, elementos, apoios, forcas)
    |       +-- elmlib() -> elmt0X()   le propriedades de material
    |
    +-- [2] numeq()    Numera as equacoes livres
    |
    +-- [3] profil()   Calcula o perfil skyline de K
    |
    +-- [4] pload()    Monta vetor global de forcas f
    |
    +-- [5] pform()    Monta K global e corrige f
    |       +-- loop nos elementos
    |       +-- elmlib() -> elmt0X()   calcula Ke local
    |       +-- addstf()               acumula Ke->K e fe->f
    |
    +-- [6] pcg()      Resolve K·u = f (Gradiente Conjugado)
    |       +-- matvec()   produto K·v no formato skyline
    |
    +-- [7] wdata()    Escreve resultados no arquivo de saida
```

### Descrição de cada sub-rotina

| Sub-rotina | Função |
|-----------|--------|
| `contr` | Controlador: lê parâmetros, gerencia memória, chama todas as etapas |
| `rdata` | Leitura completa do arquivo `.dat` |
| `numeq` | Converte flags de restrição em números de equação sequenciais |
| `profil` | Determina o perfil de esparsidade da matriz K |
| `pload` | Copia forças nodais para o vetor global |
| `pform` | Loop nos elementos: monta K e corrige f |
| `addstf` | Adiciona contribuição local (Ke, fe) à estrutura global |
| `elmlib` | Despachante: direciona para o elemento correto via `iel` |
| `elmt01..04` | Cálculo de Ke e fe para cada tipo de elemento |
| `actcol` | Solver direto LtDL (disponível, comentado no fluxo principal) |
| `pcg` | Solver iterativo: Gradiente Conjugado Precondicionado |
| `matvec` | Produto matriz-vetor K·x no formato skyline |
| `lku` | Produto s·u para forças internas locais |
| `dot` | Produto interno de dois vetores |
| `mem` | Verificação de memória disponível |
| `mzero` / `azero` | Zeragem de vetores inteiro e real |
| `wdata` | Escrita dos resultados no arquivo de saída |

### Gestão de memória em FORTRAN 77

Uma das características mais marcantes do código é o uso de **um único vetor de memória** `m(npos)`. Em FORTRAN 77 não existe alocação dinâmica (`allocate`), então toda a memória é gerenciada manualmente com ponteiros inteiros `i1, i2, i3...`:

```
Layout do vetor m/a:
|   e   |   ie  |   ix  |   id  |   x   |   f   |   u   | jdiag | am | ...
  i1      i2      i3      i4      i5      i6      i7      i8
```

O vetor é declarado como `integer m(npos)` e também como `real*8 a(1)` via `equivalence`, permitindo acessar o mesmo espaço de memória como inteiro ou real conforme necessário.

---

## 5. Tipos de elementos disponíveis

### Elemento 1 e 2 — Triângulo Linear T3

```
    3
   / \
  /   \
 1-----2
```

- **3 nós**, 2 GDL por nó → **6 GDL totais**
- Funções de forma **lineares** em x e y
- Integração **analítica** (área = det(J)/2)
- Deformações **constantes** dentro do elemento
- Requer malhas **mais refinadas** para boa precisão

**Funções de forma:**
```
N1 = 1 - x/a - y/b
N2 = x/a
N3 = y/b
```

### Elemento 3 e 4 — Quadrilátero Bilinear Q4

```
4---------3
|         |
|         |
1---------2
```

- **4 nós**, 2 GDL por nó → **8 GDL totais**
- Funções de forma **bilineares** em coordenadas isoparamétricas (r, s)
- Integração **numérica** de Gauss 2x2 (4 pontos)
- Melhor desempenho por elemento do que o T3
- Sensível à distorção: prefira elementos próximos ao quadrado

**Funções de forma isoparamétricas:**
```
N1 = (1+r)(1+s)/4    N2 = (1-r)(1+s)/4
N3 = (1-r)(1-s)/4    N4 = (1+r)(1-s)/4
```

### Transformação isoparamétrica

O Q4 é definido em coordenadas naturais (r, s) no intervalo [-1, +1], mas as integrações são feitas no espaço físico (x, y). A **matriz Jacobiana** faz esta transformação:

```
    | dx/dr   dy/dr |   | SUM(dNi/dr * xi)   SUM(dNi/dr * yi) |
J = |               | = |                                       |
    | dx/ds   dy/ds |   | SUM(dNi/ds * xi)   SUM(dNi/ds * yi) |
```

O determinante det(J) converte a área do elemento natural para o físico. Se **det(J) <= 0**, o elemento está distorcido ou com nós em ordem horária — o programa para com mensagem de erro.

---

## 6. Formato do arquivo de entrada `.dat`

### Visão geral da estrutura

```
LINHA 1:   nnode  numel  numat  nen  ndf  ndm
BLOCO A:   Materiais       (numat grupos de duas linhas)
BLOCO B:   Coordenadas     (nnode linhas)
BLOCO C:   Conectividade   (numel linhas)
BLOCO D:   Apoios          (linhas ate 0 0 0)
BLOCO E:   Forcas nodais   (linhas ate 0 0.0 0.0)
```

### Linha de cabeçalho

```
nnode  numel  numat  nen  ndf  ndm
```

| Variável | Significado | Valores típicos |
|----------|-------------|-----------------|
| `nnode` | Número total de nós | qualquer inteiro > 0 |
| `numel` | Número total de elementos | qualquer inteiro > 0 |
| `numat` | Número de materiais distintos | >= 1 |
| `nen` | Nós por elemento | `3` (T3) ou `4` (Q4) |
| `ndf` | Graus de liberdade por nó | `2` (sempre, para 2D) |
| `ndm` | Dimensão do problema | `2` (sempre) |

### Bloco A — Materiais

Para cada material, duas linhas:

```
ma  iel
E  nu  [espessura]
```

- `ma` = índice do material (1, 2, ...)
- `iel` = tipo de elemento (1 a 4 conforme tabela da seção 1)
- `E` = módulo de Young
- `nu` = coeficiente de Poisson
- `espessura` = necessária apenas para EPT (`iel = 2` ou `iel = 4`)

### Bloco B — Coordenadas nodais

Uma linha por nó:

```
no  x  y
```

> O programa foi escrito para ler três coordenadas (x, y, z), mas ignora z. Você pode incluir `0.0` ou omitir.

### Bloco C — Conectividade dos elementos

Uma linha por elemento:

```
nel  n1  n2  n3  [n4]  ma
```

- `nel` = número do elemento
- `n1..n3` (T3) ou `n1..n4` (Q4) = nós em ordem **anti-horária**
- `ma` = índice do material

> **IMPORTANTE:** A ordem dos nós deve ser **anti-horária**. Ordem horária produz det(J) <= 0 e o programa para com erro.

```
Correto (anti-horario):    Errado (horario):
    3                          1
   / \                        / \
  /   \                      /   \
 1-----2                    3-----2
```

### Bloco D — Condições de contorno

Lido em loop até encontrar `k <= 0`:

```
k  id1  id2
```

- `k` = número do nó restringido
- `id1` = flag para ux: `1` = restrito (deslocamento zero), `0` = livre
- `id2` = flag para uy: `1` = restrito (deslocamento zero), `0` = livre

Terminar com:
```
0  0  0
```

### Bloco E — Forças nodais

Lido em loop até encontrar `k <= 0`:

```
k  fx  fy
```

Terminar com:
```
0  0.0  0.0
```

---

## 7. Exemplos resolvidos

### Exemplo 1 — `Ex1.dat`: Chapa quadrada com T3 em EPD

**Problema:** Chapa 1x1, engastada no lado esquerdo, carga distribuída para baixo no lado direito.

```
  y
  |
  |  +---------------------------+  <- forcas f = -0.10 (nos 2..10)
  |  |                           |  <- forca  f = -0.05 (no 11, canto)
  |  |   121 nos, 200 elementos  |
  |  |   E = 1.0,  v = 0.30      |
  |  |   Triangulos T3, EPD      |
  |  |                           |
  |  |                           |
  |  +---------------------------+
  | apoio (nos 1,12,23,...,111)
  +--------------------------------------------> x
    ux = uy = 0
```

**Cabeçalho:**
```
121  200    1    3    2    2
```
→ 121 nós · 200 elementos · 1 material · 3 nós/elem · 2 GDL/nó · 2D

**Material (EPD, T3):**
```
  1    1
  1.0  0.30
```

**Condições de contorno (lado esquerdo engastado):**
```
  1    1    1
 12    1    1
 23    1    1
  ...
111    1    1
  0    0    0
```

**Carregamento:**
```
  2     0.0    -0.10
  3     0.0    -0.10
  ...
 10     0.0    -0.10
 11     0.0    -0.05    <- no de canto recebe metade da carga
  0     0.0     0.0
```

> **Por que o nó de canto recebe metade da carga?**
> Ao transformar uma carga distribuída uniforme em forças nodais equivalentes, um nó interno da borda recebe a contribuição de dois elementos adjacentes (carga completa), enquanto o nó de extremidade recebe de apenas um elemento (meia carga). Esta é a representação correta da carga distribuída no MEF.

---

### Exemplo 2 — `Ex2.dat`: Mesma chapa com Q4 em EPD

Mesmo domínio, mas agora com 100 elementos quadriláteros (Q4).

```
121  100    1    4    2    2
  1    3
  1.0  0.30
```

**Conectividade do elemento 1:**
```
  1    1   12   13    2    1
```
→ Elemento 1, nós 1→12→13→2 (anti-horário), material 1.

```
  No 12---------No 13
    |               |
    |   Elemento 1  |
    |               |
  No 1----------No 2
```

---

## 8. Formato do arquivo de saída

O arquivo de saída é estruturado em blocos para compatibilidade com pós-processadores:

```
coor  121
         1   0.00000E+00   1.00000E+00   0.00000E+00
         2   1.00000E-01   1.00000E+00   0.00000E+00
    ...
elem       200
         1         3    12    13     1
    ...
nosc  2
         1  -1.23456E-03  -4.56789E-03
    ...
nvec
         1  -1.23456E-03  -4.56789E-03   0.00000E+00
    ...
end
```

| Bloco | Conteúdo |
|-------|----------|
| `coor` | Coordenadas de todos os nós (x, y, z=0) |
| `elem` | Conectividade de todos os elementos |
| `nosc` | Deslocamentos nodais calculados (ux, uy) |
| `nvec` | Mesmo conteúdo de `nosc` com z=0 extra (para visualizadores 3D) |
| `end`  | Marcador de fim de arquivo |

Para nós **livres**, `nosc` contém o deslocamento calculado. Para nós **restritos**, contém o deslocamento prescrito (zero, na maioria dos casos).

---

## 9. O solver PCG explicado

O programa usa o **Gradiente Conjugado Precondicionado (PCG)** para resolver K·u = f. Este é um solver **iterativo**: ele refina uma solução aproximada a cada passo, até que o erro caia abaixo de uma tolerância.

### Por que iterativo e não direto?

O código também contém a sub-rotina `actcol` (solver direto LtDL, comentado). A tabela abaixo compara os dois:

| Característica | Solver Direto (LtDL) | Solver Iterativo (PCG) |
|----------------|----------------------|------------------------|
| Memória | Proporcional ao skyline | Apenas alguns vetores |
| Custo por solução | Depende do perfil | Depende da convergência |
| Robustez | Alta | Depende do precondicionador |
| Melhor para | Malhas pequenas/médias | Malhas grandes |

### O algoritmo PCG passo a passo

```
1. Precondicionador:   M = diag(K)   (diagonal de K)

2. Inicializar:
      x = 0                 (chute inicial)
      r = f - K*x = f       (residuo: quanto falta para equilibrar)
      d = M^{-1} * r        (direcao inicial de busca)

3. Loop ate convergir:
      z     = K * d                  (mapeia direcao pela rigidez)
      alpha = (r · d) / (d · z)      (comprimento do passo otimo)
      x     = x + alpha * d          (atualiza solucao)
      r     = r - alpha * z          (atualiza residuo)
      z     = M^{-1} * r             (aplica precondicionador)
      beta  = (r · z)_novo / (r · z)_antigo
      d     = z + beta * d           (nova direcao conjugada)

4. Parar quando  |r · z| < tol^2 * |r · z|_inicial
```

### Parâmetros no código

```fortran
tol   = 1.0d-07   ! tolerancia de convergencia
maxit = 1000      ! maximo de iteracoes
```

---

## 10. Armazenamento Skyline explicado

A matriz de rigidez K é **simétrica e esparsa**. O **formato skyline (perfil de coluna)** guarda, para cada coluna j, apenas os coeficientes desde o primeiro não-zero até a diagonal.

### Exemplo visual

```
Matriz K 5x5 (zeros omitidos):

       col1  col2  col3  col4  col5
  lin1 [K11]
  lin2 [K21  K22]
  lin3        K32  K33]
  lin4             K43  K44]
  lin5                  K54  K55]

Armazenamento em vetor (am):
  am = [K11 | K21 K22 | K32 K33 | K43 K44 | K54 K55]
         ^     ^   ^     ^   ^     ^   ^     ^   ^
        diag  col2 diag col3 diag col4 diag col5 diag

Vetor jdiag (posicao da diagonal de cada coluna):
  jdiag = [1,  3,  5,  7,  9]
```

### Por que o skyline funciona bem para MEF?

Dois coeficientes K(i,j) são **zero** se os nós i e j **nunca compartilham um elemento**. A altura do perfil depende da diferença máxima de numeração entre nós do mesmo elemento — por isso, uma boa numeração de nós reduz o tamanho do skyline.

### A sub-rotina `profil`

Percorre todos os elementos e, para cada par de equações (kk, ll) do mesmo elemento, registra a altura necessária na coluna `max(kk, ll)`:

```fortran
m = max0(kk, ll)
jdiag(m) = max0(jdiag(m), iabs(kk-ll))
```

Depois converte alturas em posições absolutas da diagonal:

```fortran
jdiag(i) = jdiag(i) + jdiag(i-1) + 1
```

---

## 11. Referências

1. **Fish, J.; Belytschko, T.** — *A First Course in Finite Elements*. Wiley, 2007.
   *(Recomendado como primeira leitura — didático e acessível)*

2. **Hughes, T.J.R.** — *The Finite Element Method*. Dover, 2000.
   *(Referência clássica, rigorosa)*

3. **Zienkiewicz, O.C.; Taylor, R.L.** — *The Finite Element Method*, 6a ed. Butterworth-Heinemann, 2005.
   *(Enciclopédia do MEF)*

4. **Bathe, K.J.** — *Finite Element Procedures*. Prentice Hall, 1996.
   *(Disponível gratuitamente no site do autor: web.mit.edu/kjb/www)*

5. **Barrett, R. et al.** — *Templates for the Solution of Linear Systems*. SIAM, 1994.
   *(Referência para o algoritmo PCG — disponível gratuitamente online)*

---

*Disciplina de Elementos Finitos · Engenharia Civil*
