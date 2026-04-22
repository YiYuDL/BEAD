


import os
from agent.agent import BEADAgent

# (Optional) Configure your API Key here
# os.environ["OPENAI_API_KEY"] = "your-api-key-here"

def BEAD():
    print("=====================================================")
    print("      Initializing BEAD Multi-Agent Orchestrator     ")
    print("=====================================================\n")

    # Instantiate the BEAD Agent
    bead_system = BEADAgent(model_name="gpt-4o")
    
    # Experimental initialization parameters (can be updated to use argparse or frontend input later)
    current_target = "PRMT5"
    core_scaffold = "c1cc(N)ccc1"
    tox_filter = "hERG liabilities"
    
    print(f"[*] Target: {current_target}")
    print(f"[*] Scaffold Constraint: {core_scaffold}")
    print(f"[*] Negative Filters: {tox_filter}\n")
    print("Type 'exit' or 'quit' to end the session.\n")

    # Start the interactive terminal conversation
    turn_counter = 1
    while True:
        try:
            user_input = input(f"\n[You - Turn {turn_counter}]: ")
            if user_input.lower() in ['exit', 'quit']:
                print("Ending BEAD session. Goodbye!")
                break
            if not user_input.strip():
                continue

            print("[BEAD is thinking...]")
            
            # Call the agent's chat method
            response = bead_system.chat(
                user_input=user_input,
                target_name=current_target,
                scaffold_smiles=core_scaffold,
                negative_filters=tox_filter
            )
            
            print(f"\n[BEAD Agent]:\n{response}\n")
            print("-" * 50)
            
            turn_counter += 1

        except KeyboardInterrupt:
            print("\nSession interrupted by user. Exiting...")
            break
        except Exception as e:
            print(f"\n[Error occurred]: {str(e)}")
            break

if __name__ == "__main__":
    BEAD()
