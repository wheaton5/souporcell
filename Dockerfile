FROM cumulusprod/souporcell:2021.03

WORKDIR /tmp
RUN rm -rf /opt/souporcell/troublet/target/release/* \
        && git clone -b check-posteriors https://github.com/johnyaku/souporcell.git . \
        && cd troublet/ \
        && cargo build \
        && cp -r /tmp/troublet/target/debug/* /opt/souporcell/troublet/target/release/ \
        && rm -rf /tmp
