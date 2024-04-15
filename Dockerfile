FROM continuumio/miniconda3
RUN conda update -n base -c defaults conda
RUN conda install -n base -c conda-forge mamba 
RUN mamba update -n base mamba 
RUN mamba create -n test_env -c bioconda python 3.11 biopython bbmap 
RUN echo "source activate test_env" > ~/.bashrc 
ENV PATH /opt/conda/envs/test_env/bin:$PATH 
RUN pip install pandas pysam matplotlib argparse 