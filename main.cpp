/* CSCI 200: Final Project - DNA Analyzer
 *
 * Author: Olivia Tallent
 * 
 * Dataset resource: https://www.kaggle.com/datasets/nageshsingh/dna-sequence-dataset/data
 *
 * Pulls DNA data from input files and create a visual analyzer to view similarities and clusters along strands
 * Press up and down arrows to navigate between strands
 * Press left and right arrows to scroll a singular strand
*/

#include "DNAStrand.h"
#include "Protein.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

vector<DNAStrand> readFile(string animalName) {
    ifstream fin("datasets/"+ animalName + ".txt");
    // check if there is an error
    if (fin.fail()) {
        cerr <<  "Error opening \'" + animalName + ".txt\' file";
        return {};
    }

    string header;
    fin >> header >> header;

    vector<DNAStrand> animal;
    // iterate through the files
    while(!fin.eof()) {
        string dnaLine;
        int classNum;
        fin >> dnaLine;
        fin >> classNum;

        // create dna strand
        DNAStrand dna = DNAStrand(animalName, dnaLine, classNum);
        animal.push_back(dna);
    }

    return animal;
}

int main() {
    int strandIndex = 0;

    // set up values for the key
    string nucleotides[4] = {"Adenosine", "Thymine", "Cytosine", "Guanine"};
    vector<sf::Color> nucleotideColors = {sf::Color::Red, sf::Color::Blue, sf::Color::Green, sf::Color::Yellow};

    string proteins[21] = {"Phenylalanine","Leucine","Serine","Tyrosine","Stop","Cysteine","Tryptophan","Proline","Histidine","Glutamine","Arginine","Isoleucine","Methionine","Threonine","Asparagine","Lysine","Valine","Alanine","Aspartate","Glutamate","Glycine"};
    vector<sf::Color> proteinColors = {sf::Color::Green, sf::Color::Blue, sf::Color::Yellow, sf::Color::Magenta, sf::Color::Red, sf::Color(35, 84, 17), sf::Color(141, 186, 224), sf::Color::Cyan, sf::Color(210, 250, 211), sf::Color(247, 233, 151), sf::Color(207, 133, 6), sf::Color(168, 63, 176), sf::Color(122, 40, 57), sf::Color(119, 71, 161), sf::Color(69, 135, 111), sf::Color(41, 0, 92), sf::Color(140, 106, 11), sf::Color(77, 16, 29), sf::Color(94, 138, 135), sf::Color::White, sf::Color(205, 255, 97)};
    
    string animal1;
    string animal2;

    // have user select animals
    do {
        cout << "Enter first animal: (chimpanzee, human, or dog) ";
        cin >> animal1;
    } while (animal1 != "chimpanzee" && animal1 != "human" && animal1 != "dog");
    do {
        cout << "Enter second animal: (chimpanzee, human, or dog) ";
        cin >> animal2;
    } while ((animal2 != "chimpanzee" && animal2 != "human" && animal2 != "dog") || animal2 == animal1);

    vector<DNAStrand> chimpanzee = readFile(animal1);
    vector<DNAStrand> dog = readFile(animal2);

    // Find size of display (min size)
    size_t end = chimpanzee.size();
    if (dog.size() < end) {
        end = dog.size();
    }

    // Create Window
    sf::Vector2u windowSize(996, 500);
    sf::RenderWindow window( sf::VideoMode( windowSize ), "DNA Analyzer" );

    int scrollPos = 0;

    while( window.isOpen() ) {
        window.clear(sf::Color(0, 0, 0));

        // pull data about comparisons from class algorithms
        vector<int> similarityClusters = chimpanzee.at(strandIndex).findClusters(dog.at(strandIndex));
        vector<int> similarityClustersP = chimpanzee.at(strandIndex).findProteinClusters(dog.at(strandIndex));
        double similiarityPercentage = chimpanzee.at(strandIndex).compareDNA(dog.at(strandIndex));
        double similiarityPercentageP = chimpanzee.at(strandIndex).compareProteins(dog.at(strandIndex));

        // display nucleotides
        chimpanzee.at(strandIndex).drawNucleotides(window, sf::Vector2f(0, 90), scrollPos);
        dog.at(strandIndex).drawNucleotides(window, sf::Vector2f(0, 115), scrollPos);
        chimpanzee.at(strandIndex).highlightNucleotideClusters(window, sf::Vector2f(10, 140), scrollPos, similarityClusters);

        // display proteins
        chimpanzee.at(strandIndex).drawProteins(window, sf::Vector2f(0, 195), scrollPos);
        dog.at(strandIndex).drawProteins(window, sf::Vector2f(0, 220), scrollPos);
        chimpanzee.at(strandIndex).highlightProteinClusters(window, sf::Vector2f(0, 245), scrollPos, similarityClustersP);

        // display text
        sf::Font myFont;
        if( !myFont.openFromFile( "datasets/arial.ttf" ) )
            return -1;
        sf::Text title( myFont );
        title.setString( "dna strand comparison: " + animal1 + " vs " + animal2 + " (strand #" + to_string(strandIndex) + ")");
        title.setPosition( sf::Vector2f(10.f, 0.f) );
        title.setFillColor( sf::Color::White );
        window.draw( title ); 

        // nucleotide header text
        sf::Text subtitle1( myFont );
        subtitle1.setString( "nucleotide clusters: ");
        subtitle1.setCharacterSize(25);
        subtitle1.setPosition( sf::Vector2f(10.f, 40.f) );
        subtitle1.setFillColor( sf::Color::White );
        window.draw( subtitle1 ); 
        sf::Text similarity1( myFont );
        similarity1.setString( "overall similarity: " + to_string(similiarityPercentage) + "%");
        similarity1.setCharacterSize(15);
        similarity1.setPosition( sf::Vector2f(10.f, 65.f) );
        similarity1.setFillColor( sf::Color::White );
        window.draw( similarity1 ); 

        // protein header text
        sf::Text subtitle2( myFont );
        subtitle2.setString( "protein clusters: ");
        subtitle2.setCharacterSize(25);
        subtitle2.setPosition( sf::Vector2f(10.f, 145.f) );
        subtitle2.setFillColor( sf::Color::White );
        window.draw( subtitle2 ); 
        sf::Text similarity2( myFont );
        similarity2.setString( "overall similarity: " + to_string(similiarityPercentageP) + "%");
        similarity2.setCharacterSize(15);
        similarity2.setPosition( sf::Vector2f(10.f, 170.f) );
        similarity2.setFillColor( sf::Color::White );
        window.draw( similarity2 ); 

        // key title
        sf::Text keytitle( myFont );
        keytitle.setString( "key:");
        keytitle.setCharacterSize(25);
        keytitle.setPosition( sf::Vector2f(10.f, 255.f) );
        keytitle.setFillColor( sf::Color::White );
        window.draw( keytitle ); 

        // create the key for nucleotides
        for (size_t i = 0; i < 4; i++) {
            sf::RectangleShape rect;
            rect.setSize(sf::Vector2f(15, 15));
            rect.setFillColor(nucleotideColors.at(i));
            rect.setPosition(sf::Vector2f(10, 300 + (float)i * 20.f));
            window.draw( rect ); 

            sf::Text keyItem( myFont );
            keyItem.setString(nucleotides[i]);
            keyItem.setCharacterSize(15);
            keyItem.setPosition( sf::Vector2f(30.f, 300 + (float)i * 20.f) );
            keyItem.setFillColor( sf::Color::White );
            window.draw( keyItem ); 
        }

        // create the key for proteins by column
        for (size_t i = 0; i < 9; i++) {
            sf::RectangleShape rect;
            rect.setSize(sf::Vector2f(15, 15));
            rect.setFillColor(proteinColors.at(i));
            rect.setPosition(sf::Vector2f(150, 300 + (float)i * 20.f));
            window.draw( rect ); 

            sf::Text keyItem( myFont );
            keyItem.setString(proteins[i]);
            keyItem.setCharacterSize(15);
            keyItem.setPosition( sf::Vector2f(170.f, 300 + (float)i * 20.f) );
            keyItem.setFillColor( sf::Color::White );
            window.draw( keyItem ); 
        }
        for (size_t i = 9; i < 18; i++) {
            sf::RectangleShape rect;
            rect.setSize(sf::Vector2f(15, 15));
            rect.setFillColor(proteinColors.at(i));
            rect.setPosition(sf::Vector2f(300, 300 + (float)(i-9) * 20.f));
            window.draw( rect ); 

            sf::Text keyItem( myFont );
            keyItem.setString(proteins[i]);
            keyItem.setCharacterSize(15);
            keyItem.setPosition( sf::Vector2f(320.f, 300 + (float)(i-9) * 20.f) );
            keyItem.setFillColor( sf::Color::White );
            window.draw( keyItem ); 
        }
        for (size_t i = 18; i < 21; i++) {
            sf::RectangleShape rect;
            rect.setSize(sf::Vector2f(15, 15));
            rect.setFillColor(proteinColors.at(i));
            rect.setPosition(sf::Vector2f(450, 300 + (float)(i-18) * 20.f));
            window.draw( rect ); 

            sf::Text keyItem( myFont );
            keyItem.setString(proteins[i]);
            keyItem.setCharacterSize(15);
            keyItem.setPosition( sf::Vector2f(470.f, 300 + (float)(i-18) * 20.f) );
            keyItem.setFillColor( sf::Color::White );
            window.draw( keyItem ); 
        }

        window.display();

        // close event
        while( const std::optional event = window.pollEvent() ) {
            if( event->is<sf::Event::Closed>() ) {
                window.close();
            }
            if (event->is<sf::Event::KeyPressed>()){
                const sf::Event::KeyPressed* keyEvent = event->getIf<sf::Event::KeyPressed>();
                // manage up/down and left/right scrolling
                if (keyEvent->code == sf::Keyboard::Key::Right) {
                    scrollPos++;
                } else if (keyEvent->code == sf::Keyboard::Key::Left && scrollPos > 0){
                    scrollPos--;
                } else if (keyEvent->code == sf::Keyboard::Key::Up && strandIndex > 0) {
                    strandIndex--;
                    scrollPos = 0;
                } else if (keyEvent->code == sf::Keyboard::Key::Down && strandIndex < (int)end - 1) {
                    strandIndex++;
                    scrollPos = 0;
                }
            }
        }
    }
    return 0;
}