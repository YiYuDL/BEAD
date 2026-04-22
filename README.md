<img width="1536" height="45" alt="image" src="https://github.com/user-attachments/assets/c7a4e1de-a022-45ca-a8a5-fa5dc2424d59" /># BEAD

Implementation of the Paper "[BEAD: a Human-in-the-Loop Agentic framework for explainable and prospective drug discovery](xxx)" (Unpublished) by Yi Yu et al.. BEAD is a novel LLM-based agent that Bridges the gap between Experiment- and Artificial intelligence-based Drug discovery.

<img src="example/Figure 1.png" width="100%" height="100%">

## Install via Anaconda
Create a new environment:
```bash
cd bead
conda env create -f environment.yml
conda activate bead
```
Install the bead package
```bash
python setup.py install
```

## Getting start
First set up your API keys in your environment.
```
export OPENAI_API_KEY=your-openai-api-key
```

## Scaffold hopping by human-in-the-loop dialogues
1st query:
```python
from bead.app import BEAD
'''Substructure Search'''
BEAD("Subsearch the structures based on Scaffold_1.mol and mol_list.csv.", "thread_003")
```
1st output is:
<img src="example/Figure 2-1.png" width="100%" height="100%">

```bash
Answer_1: The substructure search has been completed, and the results are saved in the file: `/home/user/save/search_1.csv`.
```

2nd query
```python
from bead.app import BEAD
'''Sim calculation and Rev_subsearch'''
BEAD("Filter out the structures conforming to the initial scaffold. Then calculate the similarity within the target molecule against 226192.", "thread_003")
```
1st output is:


3rd query
```python
from bead.app import BEAD
'''Sim calculation and Rev_subsearch'''
BEAD("Use two mol file (Ref.mol and New.mol) and protein file (Mat2a.pdb) to simulate the drug-target interaction.", "thread_003")
```
1st output is:
```bash
