# 🎓 Implementações dos Alunos — Elastic 2D

Esta pasta reúne as implementações individuais desenvolvidas pelos alunos como trabalho final da disciplina.

## 📁 Organização

Cada aluno deve criar **uma subpasta** com o formato `nome_sobrenome/`:

```
implementacoes_alunos/
├── joao_silva/
│   ├── README.md        ← obrigatório: descreva o que foi implementado
│   ├── src/             ← seus arquivos .f90 modificados ou novos
│   └── exemplos/        ← pelo menos um caso teste com resultado esperado
├── maria_santos/
│   └── ...
└── README.md            ← este arquivo
```

## 📝 O que deve conter o seu `README.md`

1. **Título** — nome da funcionalidade implementada
2. **Descrição** — o que foi feito em 2–4 parágrafos
3. **Formulação teórica** — equações relevantes (pode referenciar livros)
4. **Como compilar e rodar** — passo a passo
5. **Exemplo** — descrição do caso teste e resultados obtidos
6. **Verificação** — como você validou sua implementação (solução analítica, comparação com software comercial, etc.)
7. **Referências**

## 💡 Sugestões de Implementação

| Nível | Sugestão |
|-------|----------|
| ⭐ | Carregamento distribuído de superfície (pressão) |
| ⭐ | Leitura de malha no formato `.msh` do Gmsh |
| ⭐⭐ | Elemento Q4 com integração reduzida |
| ⭐⭐ | Elemento LST (Linear Strain Triangle, 6 nós) |
| ⭐⭐ | Visualizador de resultados em Python (matplotlib) |
| ⭐⭐⭐ | Elemento Q8 serendipity |
| ⭐⭐⭐ | Análise de sensibilidade de malha automatizada |
| ⭐⭐⭐ | Condições de contorno periódicas |

Converse com o professor antes de começar para alinhar o escopo.

## ⚠️ Regras

- Não modifique o código base em `../src/` — copie o que precisar para sua pasta
- Sua implementação deve ter **pelo menos um exemplo funcional**
- O `README.md` é **obrigatório** para avaliação
- Prazo: conforme cronograma da disciplina
