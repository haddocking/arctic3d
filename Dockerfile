FROM python:3.10

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /opt/software
COPY . .
RUN uv pip install . --system
WORKDIR /data
ENTRYPOINT [ "arctic3d" ]
