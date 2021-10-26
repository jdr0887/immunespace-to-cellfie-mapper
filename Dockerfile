#FROM rustlang/rust:nightly-buster as stage0
#FROM rust:1.55-buster as stage0
WORKDIR /usr/src/fuse-analysis
COPY . .
USER root
RUN id -u fa 1>/dev/null 2>&1 || (( getent group 0 1>/dev/null 2>&1 || ( type groupadd 1>/dev/null 2>&1 && groupadd -g 0 root || addgroup -g 0 -S root )) && ( type useradd 1>/dev/null 2>&1 && useradd --system --create-home --uid 1001 --gid 0 fa || adduser -S -u 1001 -G root fa ))
RUN apt-get update && apt-get install -y bash curl coreutils libc-dev libpq-dev libssl-dev openssl postgresql-client

FROM stage0 as stage1
WORKDIR /usr/src/fuse-analysis
RUN chown fa:root /usr/src/fuse-analysis
USER 1001:0
RUN cargo install --force diesel_cli --no-default-features --features postgres
RUN cargo install --force cargo-make

FROM stage1 as mainstage
WORKDIR /usr/src/fuse-analysis
USER root
EXPOSE 8000
USER 1001:0

ENTRYPOINT [ "cargo", "make", "start" ]
#ENTRYPOINT [ "sleep", "4000" ]

