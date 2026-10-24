import uuid

from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, SystemMessage
from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph import MessagesState, START, StateGraph
from langgraph.prebuilt import ToolNode, tools_condition

# Import the tools collection from the same directory and the prompt builder from the root directory
from bead.agents.make_tools import bead_toolset
from bead.prompt import build_bead_system_prompt


class BEADAgent:
    def __init__(
        self,
        model_name: str = "gpt-4.1",
        temperature: float = 0.1,
        thread_id: str | None = None,
        *,
        tools=None,
        base_url: str | None = None,
        api_key: str | None = None,
    ):
        """Initialize the BEAD agent framework."""
        model_options = {"model": model_name, "temperature": temperature}
        if base_url is not None:
            model_options["base_url"] = base_url
        if api_key is not None:
            model_options["api_key"] = api_key
        self.llm = ChatOpenAI(**model_options)
        self.tools = list(tools) if tools is not None else bead_toolset
        self.memory = MemorySaver()
        self.thread_id = thread_id or str(uuid.uuid4())
        self.config = {"configurable": {"thread_id": self.thread_id}}

    def _build_graph(self, system_prompt: str):
        """Build a non-deprecated ReAct graph for the current task context."""
        model_with_tools = self.llm.bind_tools(self.tools)

        def call_model(state: MessagesState):
            response = model_with_tools.invoke(
                [SystemMessage(content=system_prompt), *state["messages"]]
            )
            return {"messages": [response]}

        workflow = StateGraph(MessagesState)
        workflow.add_node("agent", call_model)
        workflow.add_node("tools", ToolNode(self.tools, handle_tool_errors=True))
        workflow.add_edge(START, "agent")
        workflow.add_conditional_edges("agent", tools_condition)
        workflow.add_edge("tools", "agent")
        return workflow.compile(checkpointer=self.memory)

    def run(
        self,
        user_input: str,
        target_name: str = "None specified yet",
        scaffold_smiles: str = "None specified yet",
        negative_filters: str = "None specified yet",
        recursion_limit: int = 30,
    ) -> dict:
        """Execute the ReAct loop and return the complete LangGraph state.

        Returning the state exposes AI and tool messages for reproducible agent
        evaluation while ``chat`` retains the original string-only interface.
        """

        # Dynamically assemble the system prompt for the current state
        system_prompt = build_bead_system_prompt(
            target_name=target_name,
            scaffold_smiles=scaffold_smiles,
            negative_filters=negative_filters,
        )

        # Update the graph with the latest prompt for each conversation turn.
        self.graph = self._build_graph(system_prompt)

        inputs = {"messages": [HumanMessage(content=user_input)]}
        invoke_config = {
            **self.config,
            "recursion_limit": int(recursion_limit),
        }
        return self.graph.invoke(inputs, config=invoke_config)

    def chat(
        self,
        user_input: str,
        target_name: str = "None specified yet",
        scaffold_smiles: str = "None specified yet",
        negative_filters: str = "None specified yet",
    ) -> str:
        """Execute the ReAct loop and return its final natural-language reply."""
        final_state = self.run(
            user_input=user_input,
            target_name=target_name,
            scaffold_smiles=scaffold_smiles,
            negative_filters=negative_filters,
        )
        return final_state["messages"][-1].content
