<p align="center">
  <img src="https://raw.githubusercontent.com/getbindu/create-bindu-agent/refs/heads/main/assets/light.svg" alt="bindu Logo" width="200">
</p>

<h1 align="center">BioOmni Agent</h1>
<h3 align="center">General-Purpose Biomedical AI Agent</h3>

<p align="center">
  <strong>Biomedical AI agent integrating LLM reasoning, retrieval-augmented planning, and code-based execution</strong><br/>
  Advanced tools for CRISPR screens, scRNA-seq annotation, ADMET prediction, and PDF report generation
</p>

<p align="center">
  <a href="https://github.com/Paraschamoli/biomni-agent/actions/workflows/main.yml?query=branch%3Amain">
    <img src="https://img.shields.io/github/actions/workflow/status/Paraschamoli/biomni-agent/main.yml?branch=main" alt="Build status">
  </a>
  <a href="https://img.shields.io/github/license/Paraschamoli/biomni-agent">
    <img src="https://img.shields.io/github/license/Paraschamoli/biomni-agent" alt="License">
  </a>
</p>

---

## 🎯 What is BioOmni Agent?

A comprehensive biomedical AI agent that integrates advanced LLM reasoning with specialized domain tools. Built for researchers, scientists, and healthcare professionals who need sophisticated biomedical analysis capabilities.

### Key Features
*   **🧬 CRISPR Screen Planning & Design** - Design CRISPR libraries, select optimal sgRNAs, plan validation strategies
*   **🔬 Single-Cell RNA-seq Analysis** - Annotate scRNA-seq data, identify cell types, perform differential expression
*   **💊 Drug Discovery & ADMET Prediction** - Predict absorption, distribution, metabolism, excretion, and toxicity
*   **📊 Scientific Report Generation** - Create comprehensive PDF reports with proper scientific structure
*   **🧠 Memory-Enhanced Learning** - Mem0 integration for context-aware biomedical conversations
*   **⚡ Lazy Initialization** - Fast boot times, initializes on first request

---

## 🚀 Quick Start

### Prerequisites

- Python 3.10+
- [uv](https://github.com/astral-sh/uv) package manager
- API keys for OpenRouter and Mem0 (both have free tiers)

### Installation

```bash
# Clone the repository
git clone https://github.com/Paraschamoli/biomni-agent.git
cd biomni-agent

# Create virtual environment
uv venv --python 3.12.9
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
uv sync

# Configure environment
cp .env.example .env
```

### Configuration

Edit `.env` and add your API keys:

| Key | Get It From | Required |
|-----|-------------|----------|
| `OPENROUTER_API_KEY` | [OpenRouter](https://openrouter.ai/keys) | ✅ Yes |
| `MEM0_API_KEY` | [Mem0 Dashboard](https://app.mem0.ai/dashboard/api-keys) | If you want to use Mem0 tools |

### Run the Agent

```bash
# Start the agent
uv run python -m biomni_agent

# Agent will be available at http://localhost:3773
```

### Github Setup

```bash
# Initialize git repository and commit your code
git init -b main
git add .
git commit -m "Initial commit"

# Create repository on GitHub and push (replace with your GitHub username)
gh repo create Paraschamoli/biomni-agent --public --source=. --remote=origin --push
```

---

## 💡 Usage

### Example Queries

```bash
# CRISPR screen design
"Design a CRISPR knockout screen for identifying essential genes in human cancer cell lines"

# scRNA-seq analysis
"Analyze this single-cell RNA-seq dataset to identify cell types and marker genes"

# ADMET prediction
"Predict the ADMET properties for this SMILES string: CC(C)NCC(COC1=CC=CC=C1)O"

# Experimental design
"Design an experiment to validate the role of gene X in drug resistance"
```

### Input Formats

**Plain Text:**
```
Design a CRISPR screen for identifying synthetic lethal interactions in BRCA1-deficient cells
```

**JSON:**
```json
{
  "query": "Analyze scRNA-seq data",
  "data_type": "single_cell",
  "focus": "cell_type_annotation",
  "parameters": {
    "min_cells": 10,
    "min_genes": 200
  }
}
```

### Output Structure

The agent returns structured output with:
- **Analysis Results**: Detailed biomedical analysis with interpretations
- **Experimental Plans**: Step-by-step protocols and validation strategies
- **Predictions**: ADMET properties with confidence scores and safety warnings
- **Reports**: Professional PDF documents with methods, results, and references
- **Visualizations**: Descriptive statistics and plotting recommendations

---

## 🔌 API Usage

The agent exposes a RESTful API when running. Default endpoint: `http://localhost:3773`

### Quick Start

For complete API documentation, request/response formats, and examples, visit:

📚 **[Bindu API Reference - Send Message to Agent](https://docs.getbindu.com/api-reference/all-the-tasks/send-message-to-agent)**


### Additional Resources

- 📖 [Full API Documentation](https://docs.getbindu.com/api-reference/all-the-tasks/send-message-to-agent)
- 📦 [Postman Collections](https://github.com/GetBindu/Bindu/tree/main/postman/collections)
- 🔧 [API Reference](https://docs.getbindu.com)

---

## 🎯 Skills

### biomni_agent (v1.0.0)

**Primary Capability:**
- Comprehensive biomedical research and analysis
- Integration of multiple omics data types
- Experimental design and validation planning

**Features:**
- CRISPR screen design and optimization
- Single-cell RNA-seq analysis and annotation
- ADMET property prediction from chemical structures
- PDF report generation with scientific formatting
- Memory-enhanced contextual understanding

**Best Used For:**
- Designing gene editing experiments
- Analyzing high-throughput sequencing data
- Drug discovery and safety assessment
- Generating scientific reports and publications
- Planning validation experiments

**Not Suitable For:**
- Clinical diagnosis or patient-specific medical advice
- Real-time emergency medical situations
- Regulatory compliance decisions
- Direct patient care decisions

**Performance:**
- Average processing time: ~3-8 seconds
- Max concurrent requests: 10
- Memory per request: 512MB

---

## 🐳 Docker Deployment

### Local Docker Setup

```bash
# Build and run with Docker Compose
docker-compose up --build

# Agent will be available at http://localhost:3773
```

### Docker Configuration

The agent runs on port `3773` and requires:
- `OPENROUTER_API_KEY` environment variable
- `MEM0_API_KEY` environment variable

Configure these in your `.env` file before running.

### Production Deployment

```bash
# Use production compose file
docker-compose -f docker-compose.prod.yml up -d
```

---

## 🌐 Deploy to bindus.directory

Make your agent discoverable worldwide and enable agent-to-agent collaboration.

### Setup GitHub Secrets

```bash
# Authenticate with GitHub
gh auth login

# Set deployment secrets
gh secret set BINDU_API_TOKEN --body "<your-bindu-api-key>"
gh secret set DOCKERHUB_TOKEN --body "<your-dockerhub-token>"
```

Get your keys:
- **Bindu API Key**: [bindus.directory](https://bindus.directory) dashboard
- **Docker Hub Token**: [Docker Hub Security Settings](https://hub.docker.com/settings/security)

### Deploy

```bash
# Push to trigger automatic deployment
git push origin main
```

GitHub Actions will automatically:
1. Build your agent
2. Create Docker container
3. Push to Docker Hub
4. Register on bindus.directory

---

## 🛠️ Development

### Project Structure

```
biomni-agent/
├── biomni_agent/
│   ├── skills/
│   │   └── biomni_agent/
│   │       ├── skill.yaml          # Skill configuration
│   │       └── __init__.py
│   ├── __init__.py
│   ├── __main__.py
│   ├── main.py                     # Agent entry point
│   └── agent_config.json           # Agent configuration
├── tests/
│   └── test_main.py
├── .env.example
├── docker-compose.yml
├── Dockerfile.agent
└── pyproject.toml
```

### Running Tests

```bash
make test              # Run all tests
make test-cov          # With coverage report
```

### Code Quality

```bash
make format            # Format code with ruff
make lint              # Run linters
make check             # Format + lint + test
```

### Pre-commit Hooks

```bash
# Install pre-commit hooks
uv run pre-commit install

# Run manually
uv run pre-commit run -a
```

---

## 🤝 Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Commit your changes: `git commit -m 'Add amazing feature'`
4. Push to the branch: `git push origin feature/amazing-feature`
5. Open a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Powered by Bindu

Built with the [Bindu Agent Framework](https://github.com/getbindu/bindu)

**Why Bindu?**
- 🌐 **Internet of Agents**: A2A, AP2, X402 protocols for agent collaboration
- ⚡ **Zero-config setup**: From idea to production in minutes
- 🛠️ **Production-ready**: Built-in deployment, monitoring, and scaling

**Build Your Own Agent:**
```bash
uvx cookiecutter https://github.com/getbindu/create-bindu-agent.git
```

---

## 📚 Resources

- 📖 [Full Documentation](https://Paraschamoli.github.io/biomni-agent/)
- 💻 [GitHub Repository](https://github.com/Paraschamoli/biomni-agent/)
- 🐛 [Report Issues](https://github.com/Paraschamoli/biomni-agent/issues)
- 💬 [Join Discord](https://discord.gg/3w5zuYUuwt)
- 🌐 [Agent Directory](https://bindus.directory)
- 📚 [Bindu Documentation](https://docs.getbindu.com)

---

<p align="center">
  <strong>Built with 💛 by the team from Amsterdam 🌷</strong>
</p>

<p align="center">
  <a href="https://github.com/Paraschamoli/biomni-agent">⭐ Star this repo</a> •
  <a href="https://discord.gg/3w5zuYUuwt">💬 Join Discord</a> •
  <a href="https://bindus.directory">🌐 Agent Directory</a>
</p>

#   b i o m n i - a g e n t  
 
