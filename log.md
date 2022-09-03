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
    