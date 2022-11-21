FROM continuumio/miniconda3
RUN apt update
RUN conda install -c conda-forge pymeep pymeep-extras 
RUN conda install -c conda-forge jupyterlab
WORKDIR /code/app
EXPOSE 8888
ENTRYPOINT jupyter lab --allow-root --port 8888 --ip=0.0.0.0 --no-browser
