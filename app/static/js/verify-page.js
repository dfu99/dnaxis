// Verify page canvas setup

function initVerifyCanvas(circleCoords, stapEdges, scafEdges, labels, numissues, messages) {
    var onPage = 'verify';
    var defaultRadius = 24;
    var output = fitCoordsToCanvas(circleCoords, diagram);
    var fitCoords = output[0],
        mdptx = output[1],
        mdpty = output[2];
    var circlesArr = convertCoordsToCircles(fitCoords, defaultRadius);
    translateBy(circlesArr, mdptx, mdpty);
    drawPaths(stapEdges, 15, 0.5, '#800000');
    drawPaths(scafEdges, 6, 1.0, '#000000');
    drawCircles(circlesArr, labels, numissues);

    $(document).ready(function() {
        setError(messages);
    });
}
