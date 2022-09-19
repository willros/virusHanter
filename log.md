# 2022-09-02
Starting to build the app

- Adding a lot of functionallity from COVIZ22
- Have everything stored in a SQL database?

### Sidebar
Add sidebar with possibility to choose sample to look at. 

### To plot multiple altair plots in same page:
```html 
{% block content %}

      <div id="position"></div>
      <div id="hej"></div>

      <!-- Placeholder for the tooltip -->
      <div id="vis-tooltip" class="vg-tooltip"></div>
      <!-- Render Charts -->
      <script type="text/javascript">
      function parse(url, div) {
          var opt = { mode: "vega-lite",
              renderer: "svg",
              actions: { export: true, source: false, editor: false }
              };
          vegaEmbed("#" + div, url, opt, function (error, result) {
             vegaTooltip.vegaLite(result.view, url);
             });
            }
      //Parse your Json variable here
      parse({{ histogram | safe }}, "position")
      parse({{ histogram | safe }}, "hej")

      </script>

{% endblock content%}
```

- The json object of the facet wrap plots are TOO long. It serioulsy long. Need to:
 * Save the facet wrap as png and display the png in the website instead.
 * Flask **NEEDS** to have the files inside the static folder in order to display them!
     * Moving the data to the static folder
**IT WORKS WITH THIS SETUP!!!!!**


- Run the whole app inside a docker with gunicorn and ngninx: 
    https://testdriven.io/blog/dockerizing-flask-with-postgres-gunicorn-and-nginx/
    
# 2022-09-03

Add divs and style them in the css file later. This can be useful to create a nice looking app. 
- Fix the altair plotting macro to be able to use it to plot one plot at a time. 

- Add a nice representation of the underlying tables for the data. 

- Add possibility to choose number of species to show and the cutoff threshold.

- Add so that the nextflow workflow produces the plots for contig coverage 

# 2022-09-04

- Include the read quaility and other qc stuff in the app.
     - Jinja2 expects the html templates to live within the "tempalates" folder, just like it expects images to live withing the "static" folder. Therefore, it cannot find the fastp template. Need to fix this. 
## Solved the include fastpreport problem like so:
the function copy_fastp copies the fastp report into the templates dir every time the route is called. Maybe this will not work in linux?? Will have to try. Kinda hacky solution.
    
**Include the KRONA chart**
Same way as the fastp report 

# 2022-09-07

- Add links to the contig dataframe (cat and kaiju) so that the sequence is downloaded or shown as a fasta file. 
    - Choose the name of the contig in a select window below and when submitted the user are redirected to a new url with the sequence. 
    - render_template which takes only the sequence string as input. The word-break: break-all in the css file line breaks the paragraph so the sequence is splitted. 
    
- Add possibility to choose number of species to show and the cutoff threshold. CHECK


### Compare unknown (and know) sequences between samples, when the database are big enough. 

# 2022-09-19

* Enabled a dockerized version of the app.
build the container:
```bash
# docker file:
FROM python:3.9-slim-buster
RUN mkdir -p /app
WORKDIR /app
COPY ./requirements.txt requirements.txt
RUN pip install -r requirements.txt

EXPOSE 5000

# build the image
docker build -t flask .

# run the container 
docker run -ti --rm -p 8080:5000 -v "${PWD}":/app flask flask run --host=0.0.0.0 
```
**Super important with --host=0.0.0.0**!!

ALSO WORKS WITH:
```bash
docker-compose up -d
```


