dimOtherSeries = "                  var allSeries = this.chart.series;
                                    HighchartsAddOn.mouseWasOver = true;

                                    HighchartsAddOn[this.chart.renderTo.id]=new Object();

                                    for(i=0; i<allSeries.length; i++){
                                      // Backup the series
                                      HighchartsAddOn[this.chart.renderTo.id][i] = new Object();
                                      HighchartsAddOn[this.chart.renderTo.id][i].options = owl.deepCopy(allSeries[i].options);
                                      HighchartsAddOn[this.chart.renderTo.id][i].options.color = allSeries[i].color;

                                      if(this.index != i){
                                        // Configure all other series
                                        allSeries[i].options.color = 'rgba(241,241,241,0.4)';
                                        allSeries[i].update(allSeries[i].options);
                                      }
                                    }  "

resetAllSeries = "          if(HighchartsAddOn.mouseWasOver) {
                                    var allSeries = this.chart.series;
                                    HighchartsAddOn.mouseWasOver = false;

                                    for(i=0; i<allSeries.length; i++){
                                      allSeries[i].update(HighchartsAddOn[this.chart.renderTo.id][i].options);
                                    }
                                  }  "
