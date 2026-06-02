#pragma once

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/Frame.h"
#include "GaudiKernel/MsgStream.h"

// Checks that a collection exists in a podio frame and is non-empty.
// Returns false and prints an error message if validation fails.
static bool validateCollection(
    const podio::Frame& frame,
    const std::string& collName,
    const std::string& callerName,
    MsgStream& log)
{
  try {
    const auto& coll =
        frame.get<edm4hep::SimCalorimeterHitCollection>(collName);
    log << MSG::DEBUG << "[" << callerName << "] Collection '" << collName
        << "' found with " << coll.size() << " entries." << endmsg;
    return true;
  } catch (const std::exception& e) {
    log << MSG::ERROR << "[" << callerName << "] Collection '" << collName
        << "' not found in frame: " << e.what() << endmsg;
    return false;
  }
}

// Checks that a frame parameter exists and has the expected size.
// Returns false and prints an error message if validation fails.
static bool validateParameter(
    const podio::Frame& frame,
    const std::string& paramName,
    size_t expectedSize,
    const std::string& callerName,
    MsgStream& log)
{
  auto param = frame.getParameter<std::vector<int>>(paramName);
  if (!param.has_value()) {
    log << MSG::WARNING << "[" << callerName << "] Frame parameter '"
        << paramName << "' not found. source_id will default to 0." << endmsg;
    return false;
  }
  if (param->size() != expectedSize) {
    log << MSG::ERROR << "[" << callerName << "] Frame parameter '"
        << paramName << "' has size " << param->size()
        << " but expected " << expectedSize << " (matching collection size)."
        << endmsg;
    return false;
  }
  log << MSG::DEBUG << "[" << callerName << "] Parameter '" << paramName
      << "' found with " << param->size() << " entries." << endmsg;
  return true;
}
