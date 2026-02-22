// Scaffold sequence selection logic

var use_sequences = [];
var m13_derivatives = ['m13mp18', 'p8064', 'p7308', 'p7560', '3120'];
var custom_only = ['m13mp18', 'p8064', 'p7308', 'p7560', '3120', 'custom', 'phix174'];

function update_sequences() {
    $.ajax({
        url: "/update_sequences",
        type: 'POST',
        contentType: "application/json;charset=utf-8",
        traditional: "true",
        data: JSON.stringify(use_sequences),
        dataType: "text",
        error: function(xhr, status, error) {
            setError(xhr.responseText);
        }
    })
    .done(function(result) {
    });
}

$(".button-2-Default").on('click', function () {
    if ($(this).hasClass('button-2-Chosen')) {
        $(this).removeClass('button-2-Chosen');
        $(this).addClass('button-2-Default');
        $('span#' + this.id).html('&#9674;');
        use_sequences.splice(use_sequences.indexOf(this.value), 1);

        // Re-enable m13 derivatives if needed
        if (m13_derivatives.includes(this.value)) {
            var temp_m13 = [...m13_derivatives];
            temp_m13.splice(temp_m13.indexOf(this.value), 1);
            for (let i = 0; i < temp_m13.length; i++) {
                $('button#seq-' + temp_m13[i]).removeClass('button-2-Disabled');
            }
            $('button#seq-custom').removeClass('button-2-Disabled');
        }
        // Re-enable non-custom if needed
        else if (custom_only.includes(this.value)) {
            var temp_custom = [...custom_only];
            temp_custom.splice(temp_custom.indexOf(this.value), 1);
            for (let i = 0; i < temp_custom.length; i++) {
                $('button#seq-' + temp_custom[i]).removeClass('button-2-Disabled');
            }
        }

        // Renumber
        for (let i = 0; i < use_sequences.length; i++) {
            $('span#seq-' + use_sequences[i]).html(i + 1);
        }
    } else {
        if ($(this).hasClass('button-2-Disabled')) { }
        else {
            $(this).addClass('button-2-Chosen');
            $(this).removeClass('button-2-Default');
            var lenSeq = use_sequences.length;
            $('span#' + this.id).html(lenSeq + 1);
            use_sequences.push(this.value);

            // Disable m13 derivatives
            if (m13_derivatives.includes(this.value)) {
                var temp_m13 = [...m13_derivatives];
                temp_m13.splice(temp_m13.indexOf(this.value), 1);
                for (let i = 0; i < temp_m13.length; i++) {
                    $('button#seq-' + temp_m13[i]).addClass('button-2-Disabled');
                }
                $('button#seq-custom').addClass('button-2-Disabled');
            }
            // Disable non-custom if needed
            else if (custom_only.includes(this.value)) {
                var temp_custom = [...custom_only];
                temp_custom.splice(temp_custom.indexOf(this.value), 1);
                for (let i = 0; i < temp_custom.length; i++) {
                    $('button#seq-' + temp_custom[i]).addClass('button-2-Disabled');
                }
            }
        }
    }
    update_sequences();
});

// Custom sequence textarea toggle
$('button#seq-custom').on('click', function() {
    var content = document.getElementById('custom-textarea-content');
    if (content.style.maxHeight) {
        content.style.maxHeight = null;
    } else {
        content.style.maxHeight = content.scrollHeight + "px";
    }
});

// Popup box for gear settings
$(window).on('load', function () {
    $(".trigger_popup_fricc").click(function() {
        $('.hover_bkgr_fricc').show();
    });
    $('.popupCloseButton').click(function() {
        $('.hover_bkgr_fricc').hide();
    });
});
