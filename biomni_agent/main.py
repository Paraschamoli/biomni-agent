# |---------------------------------------------------------|
# |                                                         |
# |                 Give Feedback / Get Help                |
# | https://github.com/getbindu/Bindu/issues/new/choose    |
# |                                                         |
# |---------------------------------------------------------|
#
#  Thank you users! We ❤️ you! - 🌻

"""BioOmni - A General-Purpose Biomedical AI Agent.

General-purpose biomedical AI agent integrating LLM reasoning, 
retrieval-augmented planning, and code-based execution with 
domain-specific tools for CRISPR screens, scRNA-seq annotation, 
ADMET prediction, and PDF report generation.
"""

import argparse
import asyncio
import json
import os
import sys
import traceback
from pathlib import Path
from textwrap import dedent
from typing import Any, Optional

from agno.agent import Agent
from agno.models.openrouter import OpenRouter
from agno.tools.mem0 import Mem0Tools
from bindu.penguin.bindufy import bindufy
from dotenv import load_dotenv

# Import our lightweight biomedical tools
from biomni_agent.biomedical_tools import (
    run_blast_sequence_async,
    analyze_scrnaseq_async,
    predict_admet_properties_async,
    generate_pdf_report_async,
)


# Load environment variables from .env file
load_dotenv()

# Global agent instance
agent: Agent | None = None
model_name: str | None = None
openrouter_api_key: str | None = None
mem0_api_key: str | None = None
_initialized = False
_init_lock = asyncio.Lock()


def load_config() -> dict:
    """Load agent configuration from project root."""
    # Get path to agent_config.json in project root
    config_path = Path(__file__).parent / "agent_config.json"
    
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in config file: {e}")
    except Exception as e:
        raise RuntimeError(f"Error loading config: {e}")


def get_model_name() -> str:
    """Get model name from environment or default."""
    return os.getenv("MODEL_NAME", "openai/gpt-4o")


def get_api_keys() -> tuple[Optional[str], Optional[str]]:
    """Get API keys from environment."""
    openrouter_key = os.getenv("OPENROUTER_API_KEY")
    mem0_key = os.getenv("MEM0_API_KEY")
    return openrouter_key, mem0_key


def validate_environment() -> None:
    """Validate required environment variables."""
    openrouter_key, mem0_key = get_api_keys()
    
    if not openrouter_key:
        raise ValueError(
            "OPENROUTER_API_KEY environment variable is required.\n"
            "Set it in .env file or use --api-key argument.\n"
            "Get API key from: https://openrouter.ai/keys"
        )
    
    if not mem0_key:
        print("⚠️  MEM0_API_KEY not set - memory functions will be disabled")


async def initialize_agent() -> None:
    """Initialize the agent once."""
    global agent, model_name, openrouter_api_key, mem0_api_key

    # Get settings from environment
    model_name = get_model_name()
    openrouter_api_key, mem0_api_key = get_api_keys()

    # Validate environment before initialization
    validate_environment()

    # Prepare tools list
    tools = [
        # Biomedical tools
        run_blast_sequence_async,
        analyze_scrnaseq_async,
        predict_admet_properties_async,
        generate_pdf_report_async,
    ]
    
    # Add memory tool only if API key is provided
    if mem0_api_key:
        try:
            mem0_tools = Mem0Tools(api_key=mem0_api_key)
            tools.append(mem0_tools)
        except Exception as e:
            print(f"⚠️  Failed to initialize Mem0Tools: {e}")

    # Create the biomedical agent with all our tools
    agent = Agent(
        name="BioOmni Agent",
        model=OpenRouter(
            id=model_name,
            api_key=openrouter_api_key,
            cache_response=True,
            supports_native_structured_outputs=True,
            timeout=60,
            max_retries=2,
        ),
        tools=tools,
        instructions=[dedent("""\
            You are BioOmni, a general-purpose biomedical AI agent.
            
            CORE CAPABILITIES:
            1. CRISPR Screen Planning & Design
               - Design CRISPR libraries for gene knockout/activation
               - Select optimal sgRNAs with minimal off-target effects
               - Plan control experiments and validation strategies
               - Generate screen analysis pipelines
            
            2. Single-Cell RNA-seq Analysis
               - Annotate scRNA-seq data and identify cell types
               - Perform differential expression analysis
               - Identify marker genes and cellular states
               - Generate UMAP/t-SNE visualizations (descriptive)
            
            3. Drug Discovery & ADMET Prediction
               - Predict absorption, distribution, metabolism, excretion, toxicity
               - Analyze compound properties from SMILES strings
               - Suggest lead optimization strategies
               - Identify potential safety concerns
            
            4. Experimental Hypothesis Generation
               - Generate testable biological hypotheses
               - Design experimental workflows
               - Suggest controls and validation methods
               - Predict potential outcomes
            
            5. Scientific Report Generation
               - Create comprehensive PDF reports
               - Structure findings with proper sections
               - Include methodologies and results
               - Add references and citations when possible
            
            AVAILABLE TOOLS:
            - run_blast_sequence_async: Run BLAST sequence alignment (DNA/protein)
            - analyze_scrnaseq_async: Analyze single-cell RNA-seq data
            - predict_admet_properties_async: Predict ADMET properties from SMILES
            - generate_pdf_report_async: Generate PDF reports from analysis results
            
            WORKFLOW GUIDELINES:
            1. Understand the biomedical query thoroughly
            2. Break down complex problems into executable steps
            3. Use appropriate tools for each analysis type
            4. Interpret results in biological context
            5. Generate actionable insights and next steps
            6. Create well-structured reports when requested
            
            OUTPUT FORMATS:
            - For analysis: Structured results with interpretations
            - For plans: Step-by-step experimental designs
            - For predictions: Clear probabilities with confidence levels
            - For reports: Professional PDF format with all sections
            
            SAFETY & ETHICS:
            - Always prioritize scientific accuracy
            - Flag uncertain predictions clearly
            - Consider ethical implications in biomedical research
            - Respect data privacy and confidentiality
            - Suggest validation experiments for critical findings
        """)],
        add_datetime_to_context=True,
        markdown=True,
        # show_tool_calls=True,  # Enable tool call visibility
    )
    print(f"✅ BioOmni Agent initialized with model: {model_name}")
    print(f"🔧 Loaded {len(tools)} tools")


async def run_agent(messages: list[dict[str, str]]) -> Any:
    """Run the agent with the given messages.

    Args:
        messages: List of message dicts with 'role' and 'content' keys

    Returns:
        Agent response
    """
    global agent

    if not agent:
        raise RuntimeError("Agent not initialized. Call initialize_agent() first.")
    
    try:
        # Run the agent and get response
        response = await agent.arun(messages)
        return response
    except Exception as e:
        error_msg = f"Agent execution failed: {type(e).__name__}: {e}"
        raise RuntimeError(error_msg) from e


async def handler(messages: list[dict[str, str]]) -> Any:
    """Handle incoming agent messages.

    Args:
        messages: List of message dicts with 'role' and 'content' keys
                  e.g., [{"role": "system", "content": "..."}, {"role": "user", "content": "..."}]

    Returns:
        Agent response
    """
    global _initialized

    # Validate input messages
    if not messages:
        raise ValueError("Empty messages list provided")
    
    for msg in messages:
        if not isinstance(msg, dict) or "role" not in msg or "content" not in msg:
            raise ValueError("Each message must be a dict with 'role' and 'content' keys")

    # Lazy initialization on first call
    async with _init_lock:
        if not _initialized:
            print("🔧 Initializing BioOmni Agent...")
            try:
                await initialize_agent()
                _initialized = True
            except Exception as e:
                print(f"❌ Failed to initialize agent: {e}")
                raise

    # Run the async agent
    try:
        result = await run_agent(messages)
        return result
    except Exception as e:
        print(f"❌ Agent execution error: {e}")
        raise


async def cleanup() -> None:
    """Clean up resources."""
    print("🧹 Cleaning up BioOmni resources...")
    global agent, _initialized
    agent = None
    _initialized = False


def main():
    """Run the BioOmni Agent."""
    global model_name, openrouter_api_key, mem0_api_key

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="BioOmni - Biomedical AI Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --api-key sk-xxxx --model openai/gpt-4o
  %(prog)s --api-key sk-xxxx --mem0-api-key mem0-xxxx
        """
    )
    parser.add_argument(
        "--model",
        type=str,
        default=get_model_name(),
        help="Model ID to use (default: from MODEL_NAME env var or 'openai/gpt-4o')",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=os.getenv("OPENROUTER_API_KEY"),
        help="OpenRouter API key (env: OPENROUTER_API_KEY)",
    )
    parser.add_argument(
        "--mem0-api-key",
        type=str,
        default=os.getenv("MEM0_API_KEY"),
        help="Mem0 API key (env: MEM0_API_KEY)",
    )
    parser.add_argument(
        "--config",
        type=str,
        help="Path to custom config file (optional)",
    )
    
    args = parser.parse_args()

    # Set global model name and API keys
    model_name = args.model
    openrouter_api_key = args.api_key
    mem0_api_key = args.mem0_api_key

    # Set environment variables for consistency
    if openrouter_api_key:
        os.environ["OPENROUTER_API_KEY"] = openrouter_api_key
    if mem0_api_key:
        os.environ["MEM0_API_KEY"] = mem0_api_key
    os.environ["MODEL_NAME"] = model_name

    print(f"🤖 BioOmni - Biomedical AI Agent")
    print(f"📊 Using model: {model_name}")
    print("🧬 Available tools: BLAST, scRNA-seq, ADMET, PDF Reports")
    if mem0_api_key:
        print("🧠 Mem0 memory enabled")
    else:
        print("⚠️  Mem0 memory disabled (set MEM0_API_KEY to enable)")

    try:
        # Load configuration
        if args.config:
            config_path = Path(args.config)
            if not config_path.exists():
                raise FileNotFoundError(f"Config file not found: {config_path}")
            with open(config_path, "r", encoding="utf-8") as f:
                config = json.load(f)
            print(f"📁 Using custom config: {args.config}")
        else:
            config = load_config()
        
        print(f"🌐 Server will run on: {config.get('deployment', {}).get('url', 'http://127.0.0.1:3773')}")

        # Bindufy and start the agent server
        print("🚀 Starting BioOmni server...")
        print("📝 Press Ctrl+C to stop the server")
        
        bindufy(config, handler)
        
    except KeyboardInterrupt:
        print("\n🛑 BioOmni Agent stopped by user")
    except FileNotFoundError as e:
        print(f"❌ Config error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"❌ Validation error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Cleanup on exit
        try:
            asyncio.run(cleanup())
        except Exception as e:
            print(f"⚠️  Error during cleanup: {e}")


# Bindufy and start the agent server
if __name__ == "__main__":
    main()