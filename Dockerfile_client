FROM ubuntu:latest

ENV START 0
ENV END -1
ENV THREAD 1
ENV RANGE 1500

RUN mkdir /data && \
  apt-get update && \
  apt-get install -y libx11-dev

WORKDIR /data
COPY ./god /data/god

CMD "/data/god" ${START} ${END} ${THREAD} ${RANGE}
