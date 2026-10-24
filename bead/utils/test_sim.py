from bead.utils.sim import Sim3D


REFERENCE_SMILES = "CCCCCn1c(=O)nc(NC2CCN(Cc3ccc4c(c3)OCO4)CC2)c2cc(Cl)ccc21"

QUERY_SMILES = [
    "ClC1=CC2=C(C=C1)N(C(C3=CN=CC=C3)=O)N=C2NC4CCN(CC4)CC5=CC6=C(C=C5)OCO6",
    "O=C(N1N=C(C2=C1C=CC(Cl)=C2)NC3CCN(CC4=CC(OCO5)=C5C=C4)CC3)C6CC6",
    "ClC1=CC(C(NC2CCN(CC3=CC(OCO4)=C4C=C3)CC2)=NN5C6=CC=CC=N6)=C5C=C1",
]


def main():
    for i, smiles in enumerate(QUERY_SMILES, start=1):
        similarity = Sim3D(REFERENCE_SMILES, smiles)
        print(f"molecule_{i}\t{similarity}")


if __name__ == "__main__":
    main()
