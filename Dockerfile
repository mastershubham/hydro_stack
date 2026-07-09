FROM ghcr.io/osgeo/gdal:ubuntu-small-3.12.2

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    grass \
    make \
    grass-dev \
    gdal-bin \
    libgdal-dev \
    build-essential \
    python3-pip \
    python3-dev \
    python3-venv \
    git \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /app

# Copy requirements first to leverage Docker layer caching
COPY requirements.txt .

RUN grass -c EPSG:4326 /tmp/grass_install -e && \
    grass /tmp/grass_install/PERMANENT --exec g.extension extension=r.stream.order && \
    grass /tmp/grass_install/PERMANENT --exec g.extension extension=r.stream.basins && \
    rm -rf /tmp/grass_install

RUN python3 -m venv /opt/venv --system-site-packages
RUN /opt/venv/bin/pip install --no-cache-dir -r requirements.txt

ENV PATH="/opt/venv/bin:$PATH"

COPY . .

# Create a non-root user for security
#RUN useradd -m appuser && chown -R appuser /app
#USER appuser

CMD ["python3", "hydrological_analysis.py", "--help"]
