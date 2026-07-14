# Método dos Elementos Finitos 2D – Elasticidade Ortotrópica

Este projeto implementa um código de Método dos Elementos Finitos (MEF) para análise de elasticidade bidimensional em **estado plano de tensões (EPT)**, com suporte a materiais ortotrópicos.

A versão atual estende a formulação isotrópica original para permitir propriedades direcionais distintas.

## Principais Modificações

A partir da versão anterior (isotrópica), foram introduzidas:

- Substituição do modelo constitutivo isotrópico por ortotrópico
- Inclusão de novos parâmetros de material
- Atualização das rotinas dos elementos:
  - `elmt02` (triangular T3)
  - `elmt04` (quadrilateral Q4)

## Modelo Constitutivo

Assume-se **estado plano de tensões (EPT)** com material ortotrópico.

A matriz constitutiva pode ser escrita como:

```text
D = | Ex/(1 - nuxy*nuyx)       nuxy*Ey/(1 - nuxy*nuyx)      0 |
    | nuyx*Ex/(1 - nuxy*nuyx)  Ey/(1 - nuxy*nuyx)           0 |
    | 0                        0                            Gxy |
```

com a relação de simetria:

```text
nuyx = nuxy * (Ey / Ex)
```

## Parâmetros de Entrada

Agora cada material requer:

- `Ex`: módulo de Young na direção x
- `Ey`: módulo de Young na direção y
- `nuxy`: coeficiente de Poisson
- `Gxy`: módulo de cisalhamento
- `t`: espessura

### Leitura no código

```fortran
read(nin,*) e(1), e(2), e(3), e(4), e(5)
```

## Elementos Implementados

### 1. Elemento T3 (Triangular Linear)

Subrotina: `elmt02`

- Integração analítica
- 3 nós
- Gradientes constantes
- Matriz de rigidez:

```text
Ke = B^T * D * B * (det(J)/2) * t
```

### 2. Elemento Q4 (Quadrilateral Bilinear)

Subrotina: `elmt04`

- Integração numérica Gaussiana 2x2
- 4 nós
- Jacobiana variável
- Matriz de rigidez:

```text
Ke = soma_gp [ B^T * D * B * det(J) * t ]
```

## Estrutura do Código

Arquivo principal:

- `programMEF.f`

Subrotinas relevantes:

- `elmt02` – elemento triangular T3 ortotrópico
- `elmt04` – elemento quadrilateral Q4 ortotrópico

## Observações Importantes

- A condição `det(J) > 0` é verificada (orientação anti-horária dos nós).
- A simetria material é garantida via relação entre `nuxy` e `nuyx`.
- Para recuperar o caso isotrópico:
  - `Ex = Ey = E`
  - `nuxy = nu`
  - `Gxy = E / (2 * (1 + nu))`

