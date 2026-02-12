from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from biomni_agent.main import handler


@pytest.mark.asyncio
async def test_handler_returns_response():
    """Test that handler accepts messages and returns a response."""
    messages = [{"role": "user", "content": "Hello, how are you?"}]

    # Mock the run_agent function to return a mock response
    mock_response = MagicMock()
    mock_response.run_id = "test-run-id"
    mock_response.status = "COMPLETED"

    # Mock _initialized to skip initialization and run_agent to return our mock
    with (
        patch("biomni_agent.main._initialized", True),
        patch("biomni_agent.main.run_agent", new_callable=AsyncMock, return_value=mock_response),
    ):
        result = await handler(messages)

    # Verify we get a result back
    assert result is not None
    assert result.run_id == "test-run-id"
    assert result.status == "COMPLETED"


@pytest.mark.asyncio
async def test_handler_with_multiple_messages():
    """Test that handler processes multiple messages correctly."""
    messages = [
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": "What's the weather?"},
    ]

    mock_response = MagicMock()
    mock_response.run_id = "test-run-id-2"

    with (
        patch("biomni_agent.main._initialized", True),
        patch("biomni_agent.main.run_agent", new_callable=AsyncMock, return_value=mock_response) as mock_run,
    ):
        result = await handler(messages)

    # Verify run_agent was called
    mock_run.assert_called_once_with(messages)
    assert result is not None
    assert result.run_id == "test-run-id-2"


@pytest.mark.asyncio
async def test_handler_initialization():
    """Test that handler initializes on first call."""
    messages = [{"role": "user", "content": "Test"}]

    mock_response = MagicMock()

    # Start with _initialized as False to test initialization path
    with (
        patch("biomni_agent.main._initialized", False),
        patch("biomni_agent.main.initialize_agent", new_callable=AsyncMock) as mock_init,
        patch("biomni_agent.main.run_agent", new_callable=AsyncMock, return_value=mock_response) as mock_run,
        patch("biomni_agent.main._init_lock") as mock_lock,
    ):
        # Configure the lock to work as an async context manager
        mock_lock_instance = MagicMock()
        mock_lock_instance.__aenter__ = AsyncMock(return_value=None)
        mock_lock_instance.__aexit__ = AsyncMock(return_value=None)
        mock_lock.return_value = mock_lock_instance

        result = await handler(messages)

        # Verify initialize_agent was called once
        mock_init.assert_called_once()
        # Verify run_agent was called with the messages
        mock_run.assert_called_once_with(messages)
        assert result is not None


@pytest.mark.asyncio
async def test_handler_empty_messages():
    """Test that handler rejects empty messages list."""
    with pytest.raises(ValueError, match="Empty messages list provided"):
        await handler([])


@pytest.mark.asyncio
async def test_handler_invalid_message_format():
    """Test that handler rejects invalid message format."""
    messages = [{"invalid": "format"}]

    with pytest.raises(ValueError, match="Each message must be a dict with 'role' and 'content' keys"):
        await handler(messages)
