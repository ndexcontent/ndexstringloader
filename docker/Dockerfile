# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/ndexstringloader
COPY ./ /tmp/ndexstringloader/
RUN pip install /tmp/ndexstringloader

RUN rm -rf /tmp/ndexstringloader

RUN mkdir /ndexstringloader
WORKDIR /ndexstringloader

CMD /bin/echo "To use run: ndexloadstring.py"
