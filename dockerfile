FROM continuumio/miniconda3

RUN conda update -y conda
RUN conda install -yc conda-forge python=3.10 mamba
RUN mamba install -yc conda-forge rdkit openmm openmmforcefields pdbfixer mdtraj openff-toolkit numpy biopython nglview
RUN pip install openpyxl
COPY . .
RUN pip install .
RUN python demo/run_test.py 