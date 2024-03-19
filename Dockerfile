FROM python:3.10
WORKDIR /opt/software
COPY . .
RUN pip install .
WORKDIR /data
ENTRYPOINT [ "arctic3d" ]
