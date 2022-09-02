FROM python:3.9-slim-buster

RUN mkdir -p /app
COPY . app.py /app/
WORKDIR /app
RUN pip install -r requirements.txt
EXPOSE 8080

ENTRYPOINT [ "python" ]
CMD [ "app.py" ]
