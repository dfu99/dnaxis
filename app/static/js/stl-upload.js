$(document).ready(function() {
    var previews = [];

    $('#stlup').change(function () {
        var form_data = new FormData($('#send-stl')[0]);
        $('#STLPreview').attr("src", SPINNER_URL);
        $.ajax({
            url: "/upload-stl",
            type: 'POST',
            data: form_data,
            cache: false,
            contentType: false,
            processData: false,
            dataType: "JSON",
        })
        .done(function(result) {
            previews = result;
            document.getElementById('STLPreview').src = "data:image/png;base64," + result["xview"];
        })
    });

    $('#opt_axis').change(function() {
        var e = document.getElementById("opt_axis").value;
        if (e == "X") { document.getElementById('STLPreview').src = "data:image/png;base64," + previews["xview"]; }
        if (e == "Y") { document.getElementById('STLPreview').src = "data:image/png;base64," + previews["yview"]; }
        if (e == "Z") { document.getElementById('STLPreview').src = "data:image/png;base64," + previews["zview"]; }
    });

    $('#stl_examples').change(function () {
        var val = document.getElementById('stl_examples').value;
        console.log(val);
        if (val == "none") {
            document.getElementById('STLPreview').src = "";
        } else {
            $('#STLPreview').attr("src", SPINNER_URL);
            $.ajax({
                url: "/show-example",
                type: 'POST',
                contentType: "application/json;charset=utf-8",
                traditional: "true",
                data: JSON.stringify(val),
                dataType: "JSON",
                error: function(xhr, status, error) {
                    setError(xhr.responseText);
                }
            })
            .done(function(result) {
                previews = result;
                document.getElementById('STLPreview').src = "data:image/png;base64," + result["xview"];
            })
        }
    });
});
