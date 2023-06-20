#include <TCanvas.h>
#include <TText.h>

void writeText() {
    // Create a new canvas.
    TCanvas *c1 = new TCanvas("c1", "Text Example", 200, 10, 700, 500);

    // Set margin. For example, 0.1 to the left.
    c1->SetLeftMargin(0.1);

    // Create a TText object.
    TText *t = new TText();

    // Set text align to the left (horizontal alignment = 1).
    t->SetTextAlign(11);

    // Five different lines of text.
    const int lineCount = 5;
    std::string lines[lineCount] = {
        "This is the first line",
        "This is the second line",
        "This is the third line",
        "This is the fourth line",
        "This is the fifth line"
    };

    // Loop to write the lines to the canvas.
    for(int i = 0; i < lineCount; i++) {
        // Vertical position adjusted according to line number.
        double verticalPosition = 0.8 - i * 0.1;
        t->DrawTextNDC(0.1, verticalPosition, lines[i].c_str());
    }
}

int testTLatex() {
    writeText();
    return 0;
}

