FROM python:3.10

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /opt/software
COPY . .
RUN uv sync --no-group dev --no-editable
ENV PATH="/opt/software/.venv/bin:$PATH"
WORKDIR /data
ENTRYPOINT [ "arctic3d" ]
