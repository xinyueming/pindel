#!/bin/bash
# Build script for pindel-tool Docker image

set -e

# Configuration
IMAGE_NAME="pindel-tool"
IMAGE_TAG="${1:-0.1.0}"
REGISTRY="${2:-}"

# Build
echo "Building ${IMAGE_NAME}:${IMAGE_TAG}..."
docker build -t ${IMAGE_NAME}:${IMAGE_TAG} .

# Tag for registry (if specified)
if [ -n "$REGISTRY" ]; then
    echo "Tagging for registry: ${REGISTRY}/${IMAGE_NAME}:${IMAGE_TAG}"
    docker tag ${IMAGE_NAME}:${IMAGE_TAG} ${REGISTRY}/${IMAGE_NAME}:${IMAGE_TAG}
fi

echo "Build complete."
echo ""
echo "Usage:"
echo "  docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} --help"
echo "  docker run --rm -v /path/to/data:/data ${IMAGE_NAME}:${IMAGE_TAG} pindel -f /data/ref.fa -i /data/config.txt -o /data/output"