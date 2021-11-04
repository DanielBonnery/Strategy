#' Change leaflet cluster markers to pie charts
#' 
#' @param map the map to add awesome pie chart cluster markers to
#' @param .data data for the cluster markers
#' @param .colors a vector of colors of at least the same size that nlevels(.data[[group]])
#' @param group the name of a factor variable of .data
#' @param ... further arguments to be passed to addAwesomeMarkers
#' @examples
#' data("breweries91",package="leaflet")
#' breweries91$goodbeer<-sample(as.factor(c("terrific","marvelous","culparterretaping")),nrow(breweries91),replace=T)
#' breweries91$goodbeer2<-sample(as.factor(c("terrific","marvelous","culparterretaping")),nrow(breweries91),replace=T)
#' library(leaflet)
#' library(dplyr)
#' leaflet(breweries91) %>%
#' addTiles() %>%
#' addAwesomeMarkers()
#' map<-leaflet(breweries91) %>%addTiles()
#' addpiechartclustermarkers(map,.data=breweries91,.colors=c("red","green","blue"),group="goodbeer")
#' leaflet(breweries91) %>%
#' addTiles() %>%
#' addpiechartclustermarkers(.data=breweries91,.colors=c("red","green","blue"),group="goodbeer")%>%
#' clearMarkerClusters()%>%
#' clearMarkers()%>%
#' clearGroup("goodbeer")%>%
#' addpiechartclustermarkers(.data=breweries91,
#'   .colors=c("purple","black","yellow"),
#'   group="goodbeer2")
#' map

addpiechartclustermarkers<-function(map,.data,.colors,group,...){
  
getColor <- function(x) {.colors[.data[[group]]]}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColor(1)
)

#Generate the javascript

jsscript3<-
  paste0(
"function(cluster) {
	const groups= [",paste("'",levels(.data[[group]]),"'",sep="",collapse=","),"];
	const colors= {
		groups: [",paste("'",.colors,"'",sep="",collapse=","),"],
		center:'#ddd',
		text:'black'
	};
	const markers= cluster.getAllChildMarkers();

	const proportions= groups.map(group => markers.filter(marker => marker.options.group === group).length / markers.length);
	function sum(arr, first= 0, last) {
		return arr.slice(first, last).reduce((total, curr) => total+curr, 0);
	}
	const cumulativeProportions= proportions.map((val, i, arr) => sum(arr, 0, i+1));
	cumulativeProportions.unshift(0);

	const widthgm = 2*Math.sqrt(markers.length);
	const radiusgm= 15+widthgm/2;
	const width = 2*Math.min(Math.sqrt(markers.length),5);
	const radius= 15+1.2*Math.log(markers.length)/Math.log(10)+width/2;

	const arcs= cumulativeProportions.map((prop, i) => { return {
		x   :  radius*Math.sin(2*Math.PI*prop),
		y   : -radius*Math.cos(2*Math.PI*prop),
		long: proportions[i-1] >.5 ? 1 : 0
	}});
	const paths= proportions.map((prop, i) => {
		if (prop === 0) return '';
		else if (prop === 1) return `<circle cx='0' cy='0' r='${radius}' fill='none' stroke='${colors.groups[i]}' stroke-width='${width}' stroke-alignment='center' stroke-linecap='butt' />`;
		else return `<path d='M ${arcs[i].x} ${arcs[i].y} A ${radius} ${radius} 0 ${arcs[i+1].long} 1 ${arcs[i+1].x} ${arcs[i+1].y}' fill='none' stroke='${colors.groups[i]}' stroke-width='${width}' stroke-alignment='center' stroke-linecap='butt' />`
	});

	return new L.DivIcon({
		html: `
			<svg width='60' height='60' viewBox='-30 -30 60 60' style='width: 60px; height: 60px; position: relative; top: -24px; left: -24px;' >
				<circle cx='0' cy='0' r='${radius}' stroke='none' fill='${colors.center}' />
				<text x='0' y='0' dominant-baseline='central' text-anchor='middle' fill='${colors.text}' font-size='15'>${markers.length}</text>
				${paths.join('')}
			</svg>
		`,
		className: 'marker-cluster'
	});
}")

# Generates the map.
  addAwesomeMarkers(map=map,
                    data=.data,
                    group=as.formula(paste0("~",group)),
                    icon = icons,
                    clusterOptions = markerClusterOptions(
                      iconCreateFunction =
                        JS(jsscript3)),...)}

#' Change leaflet cluster markers to pie charts
#' 
#' @param map the map to add awesome pie chart cluster markers to
#' @param .data data for the cluster markers
#' @param .colors a vector of colors of at least the same size that nlevels(.data[[group]])
#' @param group the name of a factor variable of .data
#' @param ... further arguments to be passed to addAwesomeMarkers
#' @examples
#' library(leaflet)
#' library(dplyr)
#' data("breweries91",package="leaflet")
#' breweries91$goodbeer<-sample(as.factor(c("terrific","marvelous","culparterretaping")),nrow(breweries91),replace=T)
#' leaflet(breweries91) %>%addTiles()%>%addpiechartclustermarkers2(.data=breweries91,.colors=c("red","green","blue"),group="goodbeer")


addpiechartclustermarkers2<-function(map,.data,.colors,group,clusterId="toto",...){
  
  getColor <- function(x) {.colors[.data[[group]]]}
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColor(1)
  )
  
  #Generate the javascript
  
  jsscript3<-
    paste0(
      "function(cluster) {
      const groups= [",paste("'",levels(.data[[group]]),"'",sep="",collapse=","),"];
      const colors= {
      groups: [",paste("'",.colors,"'",sep="",collapse=","),"],
      center:'#ddd',
      text:'black'
      };
      const markers= cluster.getAllChildMarkers();
      
      const proportions= groups.map(group => markers.filter(marker => marker.options.label === group).length / markers.length);
      function sum(arr, first= 0, last) {
      return arr.slice(first, last).reduce((total, curr) => total+curr, 0);
      }
      const cumulativeProportions= proportions.map((val, i, arr) => sum(arr, 0, i+1));
      cumulativeProportions.unshift(0);
      
      const widthgm = 2*Math.sqrt(markers.length);
      const radiusgm= 15+widthgm/2;
      const width = 2*Math.min(Math.sqrt(markers.length),5);
      const radius= 15+1.2*Math.log(markers.length)/Math.log(10)+width/2;
      
      const arcs= cumulativeProportions.map((prop, i) => { return {
      x   :  radius*Math.sin(2*Math.PI*prop),
      y   : -radius*Math.cos(2*Math.PI*prop),
      long: proportions[i-1] >.5 ? 1 : 0
      }});
      const paths= proportions.map((prop, i) => {
      if (prop === 0) return '';
      else if (prop === 1) return `<circle cx='0' cy='0' r='${radius}' fill='none' stroke='${colors.groups[i]}' stroke-width='${width}' stroke-alignment='center' stroke-linecap='butt' />`;
      else return `<path d='M ${arcs[i].x} ${arcs[i].y} A ${radius} ${radius} 0 ${arcs[i+1].long} 1 ${arcs[i+1].x} ${arcs[i+1].y}' fill='none' stroke='${colors.groups[i]}' stroke-width='${width}' stroke-alignment='center' stroke-linecap='butt' />`
      });
      
      return new L.DivIcon({
      html: `
      <svg width='60' height='60' viewBox='-30 -30 60 60' style='width: 60px; height: 60px; position: relative; top: -24px; left: -24px;' >
      <circle cx='0' cy='0' r='${radius}' stroke='none' fill='${colors.center}' />
      <text x='0' y='0' dominant-baseline='central' text-anchor='middle' fill='${colors.text}' font-size='15'>${markers.length}</text>
      ${paths.join('')}
      </svg>
      `,
      className: 'marker-cluster'
      });
}")

# Generates the map.
  addAwesomeMarkers(map=map,
                    data=.data,
                    label=as.formula(paste0("~",group)),
                    group=group,
                    icon = icons,
                    clusterId=rep(clusterId,nrow(.data)),
                    clusterOptions = markerClusterOptions(
                      iconCreateFunction =
                        JS(jsscript3)),...)}