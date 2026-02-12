<p align="center">
  <img src="https://raw.githubusercontent.com/getbindu/create-bindu-agent/refs/heads/main/assets/light.svg" alt="bindu Logo" width="200">
</p>

<h1 align="center">BioOmni Agent</h1>
<h3 align="center">General-Purpose Biomedical AI Agent</h3>

<p align="center">
  <strong>Biomedical AI agent integrating LLM reasoning with specialized domain tools</strong><br/>
  CRISPR screen design, scRNA-seq analysis, ADMET prediction, and scientific PDF reports
</p>

<p align="center">
  <a href="https://github.com/Paraschamoli/biomni-agent/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/Paraschamoli/biomni-agent/main.yml?branch=main" alt="Build Status">
  </a>
  <a href="https://pypi.org/project/biomni-agent/">
    <img src="https://img.shields.io/pypi/v/biomni-agent" alt="PyPI Version">
  </a>
  <img src="https://img.shields.io/badge/python-3.12+-blue.svg" alt="Python Version">
  <a href="https://github.com/Paraschamoli/biomni-agent/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/Paraschamoli/biomni-agent" alt="License">
  </a>
</p>

---

## 🎯 What is BioOmni Agent?

A comprehensive biomedical AI agent that integrates advanced LLM reasoning with specialized domain tools. Built for researchers, scientists, and bioinformaticians who need sophisticated biomedical analysis capabilities without writing code.

### Key Features
*   **🧬 CRISPR Screen Planning & Design** - Design CRISPR libraries, select optimal sgRNAs, plan validation strategies
*   **🔬 Single-Cell RNA-seq Analysis** - Annotate scRNA-seq data, identify cell types, perform differential expression
*   **💊 Drug Discovery & ADMET Prediction** - Query ChEMBL database for bioactivity data from SMILES strings
*   **🧬 BLAST Sequence Analysis** - Run NCBI BLAST searches for DNA, RNA, and protein sequences
*   **📊 Scientific Report Generation** - Create comprehensive PDF reports with proper scientific formatting
*   **🧠 Memory-Enhanced Learning** - Mem0 integration for context-aware biomedical conversations
*   **⚡ Lazy Initialization** - Fast boot times, initializes on first request
*   **🔐 Secure API Handling** - API keys from environment only, never hardcoded

---

## 🛠️ Tools & Capabilities

### Built-in Biomedical Tools
| Tool | Description | Use Case |
|------|-------------|----------|
| **`run_blast_sequence_async`** | NCBI BLAST search | Sequence homology, gene identification |
| **`analyze_scrnaseq_async`** | Single-cell RNA-seq analysis | Cell typing, clustering, marker genes |
| **`predict_admet_properties_async`** | ChEMBL database query | Drug properties, bioactivity data |
| **`generate_pdf_report_async`** | Scientific PDF generation | Research reports, documentation |

### Research Capabilities
1.  **CRISPR Screen Design** - Library design, sgRNA selection, validation strategies
2.  **Transcriptomics Analysis** - Cell annotation, differential expression, visualization
3.  **Cheminformatics** - Compound lookup, activity data, target information
4.  **Experimental Planning** - Hypothesis generation, workflow design, controls

---

> **🌐 Join the Internet of Agents**
> Register your agent at [bindus.directory](https://bindus.directory) to make it discoverable worldwide and enable agent-to-agent collaboration. It takes 2 minutes and unlocks the full potential of your agent.

---

## 🚀 Quick Start

### 1. Prerequisites

- Python 3.12+
- [uv](https://github.com/astral-sh/uv) package manager
- API keys: [OpenRouter](https://openrouter.ai/keys) (required), [Mem0](https://app.mem0.ai/dashboard/api-keys) (optional)

### 2. Clone and Setup

```bash
# Clone the repository
git clone https://github.com/Paraschamoli/biomni-agent.git
cd biomni-agent

# Set up virtual environment with uv
uv venv --python 3.12
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
uv sync
```

### 3. Configure Environment

```bash
# Copy environment template
cp .env.example .env

# Edit .env and add your API keys:
# OPENROUTER_API_KEY=sk-...  # Required for LLM
# MEM0_API_KEY=m0-...         # Optional for memory features
```

### 4. Run Locally

```bash
# Start the BioOmni agent
uv run python -m biomni_agent

# Or directly
uv run biomni-agent

# Access at: http://localhost:3773
```

### 5. Test with Docker

```bash
# Build and run with Docker Compose
docker-compose up --build

# Access at: http://localhost:3773
```

---

## 🔧 Configuration

### Environment Variables
Create a `.env` file:

```env
# Required
OPENROUTER_API_KEY=sk-...     # OpenRouter API key (get from openrouter.ai)

# Optional
MEM0_API_KEY=m0-...           # Mem0 API key for memory features
MODEL_NAME=openai/gpt-4o      # Model to use (default: openai/gpt-4o)
DEBUG=true                    # Enable debug logging
```

### Port Configuration
Default port: `3773` (can be changed in `agent_config.json`)

---

## 💡 Usage Examples

### Via HTTP API

```bash
curl -X POST http://localhost:3773/chat \
  -H "Content-Type: application/json" \
  -d '{
    "messages": [
      {
        "role": "user",
        "content": "Design a CRISPR knockout screen for identifying essential genes in human cancer cell lines. Include sgRNA selection criteria and validation strategies."
      }
    ]
  }'
```

### Sample Research Queries

```text
# CRISPR screen design
"Design a CRISPR-Cas9 screen to identify genes involved in drug resistance in melanoma"

# scRNA-seq analysis
"Analyze this single-cell RNA-seq dataset to identify cell types and marker genes: /path/to/data.h5ad"

# ADMET prediction
"Query ChEMBL for bioactivity data on this SMILES: CC(C)NCC(COC1=CC=CC=C1)O"

# BLAST search
"Run BLAST on this DNA sequence: AGCTAGCTAGCTAGCTAGCTAGCT"

# PDF report generation
"Generate a comprehensive PDF report on CRISPR-Cas9 applications in immunotherapy"
```

### Expected Output Format

```json
{
  "success": true,
  "n_cells": 12500,
  "n_genes": 18000,
  "n_clusters": 8,
  "cluster_sizes": {
    "0": 3200,
    "1": 2800,
    "2": 2100,
    "3": 1500,
    "4": 1200,
    "5": 900,
    "6": 500,
    "7": 300
  },
  "summary": "Identified 8 distinct cell clusters including T cells, monocytes, and B cells."
}
```

---

## 🐳 Docker Deployment

### Quick Docker Setup

```bash
# Build the image
docker build -f Dockerfile.agent -t biomni-agent .

# Run container
docker run -d \
  -p 3773:3773 \
  -e OPENROUTER_API_KEY=your_key_here \
  -e MEM0_API_KEY=your_mem0_key \
  --name biomni-agent \
  biomni-agent

# Check logs
docker logs -f biomni-agent
```

### Docker Compose (Recommended)

**docker-compose.yml**
```yaml
version: '3.8'
services:
  biomni-agent:
    build:
      context: .
      dockerfile: Dockerfile.agent
    ports:
      - "3773:3773"
    environment:
      - OPENROUTER_API_KEY=${OPENROUTER_API_KEY}
      - MEM0_API_KEY=${MEM0_API_KEY}
      - MODEL_NAME=${MODEL_NAME:-openai/gpt-4o}
    env_file:
      - .env
    restart: unless-stopped
    volumes:
      - ./data:/app/data
```

**Run with Compose:**
```bash
# Start with compose
docker-compose up -d

# View logs
docker-compose logs -f

# Stop
docker-compose down
```

---

## 📁 Project Structure

```
biomni-agent/
├── biomni_agent/
│   ├── __init__.py              # Package initialization
│   ├── __main__.py              # Entry point
│   ├── __version__.py           # Version information
│   ├── main.py                  # Agent implementation
│   ├── biomedical_tools.py      # Biomedical tool functions
│   ├── agent_config.json        # Bindu agent configuration
│   └── skills/
│       └── biomni/
│           ├── skill.yaml       # Skill metadata
│           └── __init__.py
├── tests/
│   ├── __init__.py
│   └── test_main.py            # Unit tests
├── pyproject.toml              # Python dependencies
├── Dockerfile.agent            # Docker build configuration
├── docker-compose.yml          # Docker Compose setup
├── README.md                   # This documentation
├── .env.example                # Environment template
└── LICENSE                     # MIT License
```

---

## 🔌 API Reference

### Health Check

```bash
GET http://localhost:3773/health
```

**Response:**
```json
{"status": "healthy", "agent": "BioOmni Agent"}
```

### Chat Endpoint

```bash
POST http://localhost:3773/chat
Content-Type: application/json

{
  "messages": [
    {"role": "user", "content": "Your biomedical research query here"}
  ]
}
```

**Response:**
```json
{
  "run_id": "abc-123",
  "content": "## CRISPR Screen Design Results...",
  "status": "COMPLETED"
}
```

For complete API documentation, visit: [Bindu API Reference](https://docs.getbindu.com/api-reference/all-the-tasks/send-message-to-agent)

---

## 🧪 Testing

### Local Testing

```bash
# Install test dependencies
uv sync --group test

# Run all tests
uv run pytest tests/

# Run with coverage
uv run pytest --cov=biomni_agent tests/

# Run specific test
uv run pytest tests/test_main.py -v
```

### Integration Test

```bash
# Start agent in background
uv run biomni-agent &

# Test API endpoint
curl -X POST http://localhost:3773/chat \
  -H "Content-Type: application/json" \
  -d '{"messages": [{"role": "user", "content": "Design a CRISPR screen for essential genes"}]}'

# Kill agent
kill %1
```

---

## 🚨 Troubleshooting

### Common Issues & Solutions

| Issue | Solution |
|-------|----------|
| "OPENROUTER_API_KEY not set" | Create `.env` file with your API key or pass via `--api-key` |
| "ModuleNotFoundError" | Run `uv sync --force` to reinstall dependencies |
| "Port 3773 already in use" | Change port in `agent_config.json` or kill process: `lsof -ti:3773 | xargs kill -9` |
| "ChEMBL client unavailable" | Install chembl-webresource-client: `uv add chembl-webresource-client` |
| "BLAST search timed out" | Use shorter sequence (<1000bp) or try again later |
| "File not found" | Ensure data path is absolute or relative to working directory |
| Docker build fails | `docker system prune -a && docker-compose build --no-cache` |
| Memory issues | Increase memory limit in Docker or reduce dataset size |

---

## 📊 Dependencies

### Core Packages

| Package | Version | Purpose |
|---------|---------|---------|
| bindu | >=2026.6.6 | Agent deployment framework |
| agno | >=2.2.0 | AI agent framework |
| openai | >=2.11.0 | LLM client |
| biopython | >=1.86 | BLAST, sequence analysis |
| scanpy | >=1.12 | Single-cell RNA-seq analysis |
| chembl-webresource-client | >=0.10.9 | ChEMBL database queries |
| reportlab | >=4.4.9 | PDF generation |
| pandas | >=2.0.0 | Data manipulation |
| python-dotenv | >=1.0.0 | Environment management |

### Development Packages

| Package | Purpose |
|---------|---------|
| pytest | Testing framework |
| ruff | Code formatting/linting |
| pre-commit | Git hooks |
| mypy | Type checking |

---

## 🤝 Contributing

We welcome contributions from the biomedical research and AI communities!

### Contribution Steps

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/amazing-feature`
3. **Make your changes** following the code style
4. **Add tests** for new functionality
5. **Run tests**: `uv run pytest`
6. **Run linters**: `uv run ruff check . && uv run ruff format .`
7. **Commit with descriptive messages**: `git commit -m 'Add amazing feature'`
8. **Push to your fork**: `git push origin feature/amazing-feature`
9. **Open a Pull Request**

### Code Style Guidelines

- Follow PEP 8 conventions
- Use type hints for all function parameters and returns
- Add docstrings for all public functions (Google style)
- Keep functions focused and small (<50 lines)
- Use async/await for I/O-bound operations
- Handle errors gracefully with descriptive messages

---

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

Copyright (c) 2026 Paras Chamoli

---

## 🙏 Credits & Acknowledgments

- **Developer:** Paras Chamoli
- **Framework:** [Bindu](https://bindus.directory) - Agent deployment platform for the Internet of Agents
- **Agent Framework:** [Agno](https://github.com/agno-agi/agno) - Lightweight AI agent toolkit
- **Bioinformatics:** [Biopython](https://biopython.org/) - NCBI BLAST integration
- **Single-Cell:** [Scanpy](https://scanpy.readthedocs.io/) - scRNA-seq analysis
- **Cheminformatics:** [ChEMBL](https://www.ebi.ac.uk/chembl/) - Drug discovery database
- **PDF Generation:** [ReportLab](https://www.reportlab.com/) - Scientific report creation

### 🔗 Useful Links

- 🌐 **Bindu Directory:** [bindus.directory](https://bindus.directory)
- 📚 **Bindu Docs:** [docs.getbindu.com](https://docs.getbindu.com)
- 🐙 **GitHub:** [github.com/Paraschamoli/biomni-agent](https://github.com/Paraschamoli/biomni-agent)
- 📖 **Documentation:** [Paraschamoli.github.io/biomni-agent](https://Paraschamoli.github.io/biomni-agent)
- 💬 **Discord:** [Bindu Community](https://discord.gg/3w5zuYUuwt)
- 🐛 **Issue Tracker:** [github.com/Paraschamoli/biomni-agent/issues](https://github.com/Paraschamoli/biomni-agent/issues)

---

## 🧬 Citation

If you use BioOmni Agent in your research, please cite:

```bibtex
@software{chamoli2026biomni,
  author = {Chamoli, Paras},
  title = {BioOmni Agent: General-Purpose Biomedical AI Agent},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/Paraschamoli/biomni-agent}
}
```

---

<p align="center">
  <strong>Built with 💛 for the biomedical research community</strong><br/>
  <em>Democratizing access to advanced bioinformatics analysis</em>
</p>

<p align="center">
  <a href="https://github.com/Paraschamoli/biomni-agent/stargazers">⭐ Star on GitHub</a> •
  <a href="https://bindus.directory">🌐 Register on Bindu</a> •
  <a href="https://github.com/Paraschamoli/biomni-agent/issues">🐛 Report Issues</a> •
  <a href="https://discord.gg/3w5zuYUuwt">💬 Join Discord</a>
</p>

> **Note:** This agent follows the Bindu pattern with lazy initialization and secure API key handling. It boots without API keys and only fails at runtime if keys are needed but not provided. BioOmni is designed for research assistance and is not a substitute for professional medical advice, clinical diagnosis, or regulatory decision-making.
