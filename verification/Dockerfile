FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential git python3 ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /evidence
CMD ["bash", "verification/clean_room_quick_verify.sh"]
