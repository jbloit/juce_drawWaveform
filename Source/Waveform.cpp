/*
  ==============================================================================

    Waveform.cpp
    Created: 3 May 2022 6:18:30am
    Author:  Julien Bloit

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Waveform.h"

//==============================================================================
Waveform::Waveform()
{
    afm.registerBasicFormats();
    thumbnail.reset( new AudioThumbnailBars(512, afm , thumbnailCache) );

    thumbnail->addChangeListener(this);

    audioFile = juce::File (juce::File::getSpecialLocation(juce::File::SpecialLocationType::tempDirectory).getChildFile("temp.ogg"));
    audioFile.replaceWithData(BinaryData::countdown_ogg, BinaryData::countdown_oggSize);

    thumbnail->setSource(new juce::FileInputSource(audioFile));

}

Waveform::~Waveform()
{
    audioFile.deleteFile();
}

void Waveform::changeListenerCallback (juce::ChangeBroadcaster* source)
{
    if (source == thumbnail.get())
    {
        repaint();
    }
}

void Waveform::paint (juce::Graphics& g)
{
    /* This demo code just fills the component's background and
       draws some placeholder text to get you started.

       You should replace everything in this method with your own
       drawing code..
    */

    g.fillAll (juce::Colours::black);   // clear the background

    g.setColour (juce::Colours::black);
    g.drawRect (getLocalBounds(), 1);   // draw an outline around the component

    g.setColour (juce::Colours::orangered);

    thumbnail->drawChannels (g, getLocalBounds(),0.0, thumbnail->getTotalLength(), 1.0f);
}

void Waveform::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..

}
