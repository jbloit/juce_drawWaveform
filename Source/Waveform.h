/*
  ==============================================================================

    Waveform.h
    Created: 3 May 2022 6:18:30am
    Author:  Julien Bloit

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class Waveform  : public juce::Component, public juce::ChangeListener
{
public:
    Waveform();
    ~Waveform() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:

    juce::AudioFormatManager afm;
    juce::AudioThumbnailCache thumbnailCache  { 5 };
    std::unique_ptr<juce::AudioThumbnail> thumbnail;
    juce::File audioFile;
    void changeListenerCallback (juce::ChangeBroadcaster* source) override;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Waveform)


};
