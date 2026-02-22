// Pathway page AJAX submit and canvas setup

$(document).ready(function() {
    $('#submit_pathway').on('click', function() {
        $.ajax({
            url: "/upload_pathway",
            type: 'POST',
            contentType: "application/json;charset=utf-8",
            traditional: "true",
            data: JSON.stringify(connections),
            dataType: "text",
            error: function(xhr, status, error) {
                setError(xhr.responseText);
            }
        })
        .done(function(result) {
            window.location.href = '/pre-submit';
        });
    });
});

function initPathwayCanvas(circleCoords, edgesData) {
    var onPage = 'pathway';
    var defaultRadius = 24;
    var output = fitCoordsToCanvas(circleCoords, diagram);
    var fitCoords = output[0],
        mdptx = output[1],
        mdpty = output[2];
    var circlesArr = convertCoordsToCircles(fitCoords, defaultRadius);
    translateBy(circlesArr, mdptx, mdpty);
    var edges = edgesData;
    drawPaths(edges, 15, 0.5, '#800000');
    drawCircles(circlesArr);

    $('#helpbtn').hover(
        function() { $('.helpimg').show(); },
        function() { $('.helpimg').hide(); }
    );
}
