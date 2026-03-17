# 📖 Guia de Git para Iniciantes

> Um guia prático para alunos de Engenharia Civil que estão usando Git pela primeira vez.

---

## O que é Git e por que usar?

Imagine que você está escrevendo sua monografia. Você salva versões como `monografia_final.docx`, `monografia_final_v2.docx`, `monografia_agora_vai.docx`... Git resolve exatamente isso — de forma elegante.

**Git** é um sistema de controle de versão. Ele:
- Guarda o histórico completo de todas as mudanças no seu código
- Permite voltar a qualquer versão anterior
- Facilita o trabalho em equipe sem "sobrescrever" o trabalho do outro
- É a ferramenta padrão da indústria de software (e cada vez mais na engenharia)

**GitHub** é um site que hospeda repositórios Git na nuvem, como se fosse um "Google Drive para código", mas muito mais poderoso.

---

## 🛠️ Instalação

### Windows
Baixe e instale o [Git for Windows](https://git-scm.com/download/win).
Durante a instalação, deixe todas as opções no padrão.

### macOS
Abra o Terminal e execute:
```bash
xcode-select --install
```
Ou instale via [Homebrew](https://brew.sh/): `brew install git`

### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install git
```

### Verificando a instalação
```bash
git --version
# Deve retornar algo como: git version 2.x.x
```

---

## ⚙️ Configuração Inicial (faça uma vez só)

```bash
git config --global user.name "Seu Nome Completo"
git config --global user.email "seu.email@exemplo.com"
```

> Use o mesmo e-mail da sua conta no GitHub!

---

## 🗺️ Conceitos Essenciais

Antes de ver os comandos, entenda os conceitos:

| Conceito | O que é |
|----------|---------|
| **Repositório (repo)** | A pasta do projeto com todo o histórico de versões |
| **Commit** | Um "snapshot" do projeto em um momento — como um ponto de save |
| **Branch** | Uma linha de desenvolvimento paralela (você trabalha sem afetar o código principal) |
| **Clone** | Baixar uma cópia completa de um repositório remoto |
| **Push** | Enviar seus commits locais para o GitHub |
| **Pull** | Baixar atualizações do GitHub para sua máquina |
| **Fork** | Criar uma cópia pessoal de um repositório de outra pessoa |
| **Pull Request (PR)** | Pedir para que suas mudanças sejam incorporadas ao repositório original |

---

## 📋 Fluxo de Trabalho Básico

Este é o ciclo que você repetirá dezenas de vezes:

```
Editar arquivos → Adicionar ao stage → Commitar → Enviar (push)
```

### Passo a passo

**1. Clone o repositório (só na primeira vez)**
```bash
git clone https://github.com/anabeatrizgonzaga/FEM-course.git
cd FEM-course
```

**2. Veja o estado atual do seu repositório**
```bash
git status
```

**3. Após editar arquivos, adicione-os ao "stage"**
```bash
# Adicionar um arquivo específico
git add nome_do_arquivo.f90

# Adicionar todos os arquivos modificados
git add .
```

**4. Crie um commit com uma mensagem descritiva**
```bash
git commit -m "Adiciona elemento Q4 com integração de Gauss 2x2"
```
> ⚠️ Escreva mensagens que descrevam **o que foi feito**, não "mudanças" ou "atualização".

**5. Envie para o GitHub**
```bash
git push origin main
```

**6. Antes de trabalhar, sempre atualize seu repositório local**
```bash
git pull origin main
```

---

## 🌿 Trabalhando com Branches

Para o trabalho final, você trabalhará em uma branch separada:

```bash
# Crie e mude para uma nova branch
git checkout -b feat/joao-silva

# Veja em qual branch você está
git branch

# Volte para a branch principal
git checkout main

# Mescle sua branch na principal (após revisão do professor)
git merge feat/joao-silva
```

**Por que branches?** Você pode desenvolver sua implementação sem "quebrar" o código principal. O professor pode revisar antes de incorporar.

---

## 🔁 Fork + Pull Request (fluxo para o trabalho)

Este é o fluxo que usaremos para as implementações individuais:

### 1. Faça um Fork
No GitHub, acesse o repositório da disciplina e clique em **Fork** (canto superior direito). Isso cria uma cópia na sua conta.

### 2. Clone o *seu* fork
```bash
git clone https://github.com/anabeatrizgonzaga/FEM-course.git
cd FEM-course
```

### 3. Adicione o repositório original como "upstream"
```bash
git remote add upstream https://github.com/anabeatrizgonzaga/FEM-course.git
```

### 4. Crie sua branch de trabalho
```bash
git checkout -b feat/joao-silva
```

### 5. Desenvolva, commite e envie
```bash
git add .
git commit -m "Implementa elemento Q8 serendipity"
git push origin feat/joao-silva
```

### 6. Abra um Pull Request
No GitHub, acesse seu fork e clique em **"Compare & pull request"**. Descreva o que você implementou.

### 7. Sincronize com o upstream (para não ficar desatualizado)
```bash
git fetch upstream
git checkout main
git merge upstream/main
```

---

## 🕵️ Inspecionando o Histórico

```bash
# Ver histórico de commits
git log

# Versão compacta e visual
git log --oneline --graph --all

# Ver o que mudou em um commit específico
git show abc1234

# Ver diferenças ainda não commitadas
git diff
```

---

## ⏪ Desfazendo Coisas

```bash
# Desfazer mudanças em um arquivo (antes do commit)
git checkout -- nome_do_arquivo.f90

# Remover arquivo do stage (sem desfazer edições)
git reset HEAD nome_do_arquivo.f90

# Voltar para o último commit (CUIDADO: descarta tudo não commitado)
git reset --hard HEAD
```

> 💡 **Regra de ouro:** Se você fez um commit, seus dados estão seguros. Git raramente perde dados commitados.

---

## ✅ Boas Práticas

### Mensagens de commit

| ✅ Bom | ❌ Ruim |
|--------|---------|
| `Corrige montagem da matriz de rigidez para Q4` | `fix` |
| `Adiciona leitura de malha no formato Gmsh` | `mudanças` |
| `Refatora módulo de integração de Gauss` | `atualização 2` |

### Outros hábitos importantes

- **Commite com frequência** — pequenas unidades de trabalho são mais fáceis de revisar e desfazer
- **Sempre faça `git pull` antes de começar a trabalhar** — evita conflitos
- **Nunca commite arquivos compilados** (`.o`, `.exe`, binários) — use um `.gitignore`
- **Uma coisa por commit** — não misture "adicionei elemento novo" com "corrigi bug no solver"

---

## 📄 Arquivo `.gitignore`

Crie um arquivo `.gitignore` na raiz do projeto para evitar commitar arquivos desnecessários:

```gitignore
# Arquivos compilados Fortran
*.o
*.mod
*.exe

# Binários
truss2d
elastic2d

# Arquivos de saída
*.dat.out
resultados/

# Arquivos de sistema
.DS_Store
Thumbs.db

# Editor
*.swp
.vscode/
```

---

## 🆘 Problemas Comuns

### "Permission denied (publickey)"
Você precisa configurar uma chave SSH. Siga [este guia do GitHub](https://docs.github.com/pt/authentication/connecting-to-github-with-ssh).

### "Merge conflict"
Acontece quando dois commits modificam a mesma linha. Git marcará os conflitos no arquivo:
```
<<<<<<< HEAD
  sua versão
=======
  versão do colega
>>>>>>> feat/outro-aluno
```
Edite o arquivo, escolha (ou combine) as versões, e depois:
```bash
git add arquivo_conflitante.f90
git commit -m "Resolve conflito de merge em element.f90"
```

### "Detached HEAD"
Você está vendo um commit antigo. Volte para a branch:
```bash
git checkout main
```

---

## 📚 Para Aprender Mais

- [Git — Documentação Oficial (PT)](https://git-scm.com/book/pt-br/v2)
- [GitHub Skills](https://skills.github.com/) — cursos interativos gratuitos
- [Learn Git Branching](https://learngitbranching.js.org/?locale=pt_BR) — visual e interativo, muito recomendado!
- [Oh Shit, Git!](https://ohshitgit.com/pt_BR) — soluções para situações de pânico

---

<p align="center">
  Dúvidas? Abra uma <strong>Issue</strong> no repositório ou fale com o professor em aula.
</p>
