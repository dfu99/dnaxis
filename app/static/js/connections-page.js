// Connections page AJAX submit and canvas setup

$(document).ready(function() {
    $('#submit_connections').on('click', function() {
        $.ajax({
            url: "/upload_connections",
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
            window.location.href = '/pathway';
        });
    });
});

function initConnectionsCanvas(circleCoords, existingConnections) {
    var onPage = 'connections';
    var defaultRadius = 24;
    var output = fitCoordsToCanvas(circleCoords, diagram);
    var fitCoords = output[0],
        mdptx = output[1],
        mdpty = output[2];
    var circlesArr = convertCoordsToCircles(fitCoords, defaultRadius);
    translateBy(circlesArr, mdptx, mdpty);
    var edges = [];
    drawCircles(circlesArr);
    applyConnections(existingConnections);

    $('#helpbtn').hover(
        function() { $('.helpimg').show(); },
        function() { $('.helpimg').hide(); }
    );
}
