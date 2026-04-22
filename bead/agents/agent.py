


import uuid
from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage
from langgraph.prebuilt import create_react_agent
from langgraph.checkpoint.memory import MemorySaver

# Import the tools collection from the same directory and the prompt builder from the root directory
from agent.make_tools import bead_toolset
from prompt import build_bead_system_prompt

class BEADAgent:
    def __init__(self, model_name: str = "gpt-4o", temperature: float = 0.1):
        """Initialize the BEAD agent framework."""
        self.llm = ChatOpenAI(model=model_name, temperature=temperature)
        self.tools = bead_toolset
        self.memory = MemorySaver()
        self.thread_id = str(uuid.uuid4())
        self.config = {"configurable": {"thread_id": self.thread_id}}

    def chat(self, user_input: str, target_name: str, scaffold_smiles: str, negative_filters: str) -> str:
        """Execute the ReAct loop for tool invocation with memory."""
        
        # Dynamically assemble the system prompt for the current state
        system_prompt = build_bead_system_prompt(
            target_name=target_name,
            scaffold_smiles=scaffold_smiles,
            negative_filters=negative_filters
        )

        # Update the graph with the latest prompt for each conversation turn
        self.graph = create_react_agent(
            self.llm, 
            tools=self.tools, 
            checkpointer=self.memory,
            state_modifier=system_prompt
        )

        inputs = {"messages": [HumanMessage(content=user_input)]}
        final_state = self.graph.invoke(inputs, config=self.config)
        
        return final_state["messages"][-1].content
