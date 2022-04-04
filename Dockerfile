FROM heat:latest

RUN mkdir -p /code/UKCP18dir /code/PreProcessedData
COPY ./launch_heat.sh /code/launch_heat.sh
RUN chmod 777 /code/launch_heat.sh
WORKDIR /code

# This is the command that will run your model
ENTRYPOINT ["/code/launch_heat.sh"]


